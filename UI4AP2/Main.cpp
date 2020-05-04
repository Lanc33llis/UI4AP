#include <SFML/Graphics.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <Windows.h>
#include <iostream>
#include <functional>
#include "resource.h"
#include "AutoPilot.h"
#include "RoundedRectangle.h"
//mainCRTStartup, CONSOLE

#define within(x, y) getGlobalBounds().contains(x, y)

sf::Image readCVGraph(std::string fileadress, sf::Text &timeNumber)
{
    cv::Mat graph = cv::Mat::zeros(cv::Size(800, 400), CV_8UC4);
    cv::Mat lines = cv::Mat::zeros(cv::Size(800, 400), CV_8UC4);

    AP::Path myPath = { AP::Waypoint{3.75, 5.5, 180}, AP::Waypoint{3, 5.5, 180}, /*AP::Waypoint{6, 4, 80}, AP::Waypoint{6, 4, 100},*/ AP::Waypoint{6.5, 7.1, 0}, AP::Waypoint{8, 7.1, 0} };
    AP::Trajectory myTrajectory = AP::TrajectoryGeneration(myPath, 2);

    double xRatio, yRatio, xLength = 16.5, yLength = 8.23;
    xRatio = graph.cols / xLength; yRatio = graph.rows / yLength;

    sf::Image buffer;

    for (std::size_t i = 0; i <= 16; i++)
    {
        cv::line(lines, cv::Point(i * xRatio, lines.rows), cv::Point(i * xRatio, 0), cv::Scalar(1, 1, 1, 255), 1, 8);
    }

    for (std::size_t i = 0; i <= 8; i++)
    {
        cv::line(lines, cv::Point(0, i * yRatio), cv::Point(lines.cols, i * yRatio), cv::Scalar(1, 1, 1, 255), 1, 8);
    }

    for (std::size_t i = 0; i < myTrajectory.size(); i++)
    {
        double c0 = myTrajectory[i].Function.Ax, c1 = myTrajectory[i].Function.Bx, c2 = myTrajectory[i].Function.Cx, c3 = myTrajectory[i].Function.Dx;

        std::function<double(double)> cubicFunction = [c0, c1, c2, c3](double x)
        {
            return ((c0 * pow(x, 3)) + (c1 * pow(x, 2)) + (c2 * x) + c3);
        };

        std::vector<cv::Point> points;

        if ((myTrajectory[i].Function.PointTwo.X - myTrajectory[i].Function.PointOne.X) < 0)
        {
            for (double g = myTrajectory[i].Function.PointOne.X; g >= myTrajectory[i].Function.PointTwo.X; g -= .005)
            {
                points.push_back(cv::Point(g * xRatio, cubicFunction(g) * yRatio));
            }
        }

        else
        {
            for (double g = myTrajectory[i].Function.PointOne.X; g <= myTrajectory[i].Function.PointTwo.X; g += .005)
            {
                points.push_back(cv::Point(g * xRatio, cubicFunction(g) * yRatio));
            }
        }

        cv::polylines(graph, points, false, cv::Scalar(255, 255, 0, 1), 3, 8, 0);
    }

    cv::Mat rest(graph.rows, graph.cols, CV_8UC3);
    cv::flip(graph, graph, 0);
    cv::flip(lines, lines, 0);

    if (fileadress.empty())
    {
        rest.setTo(cv::Scalar(240, 240, 240, 255));
        rest.setTo(cv::Scalar(0, 0, 0, 254), graph);
        rest.setTo(cv::Scalar(0, 0, 0, 255), lines);
    }

    else
    {
        rest = cv::imread(fileadress, cv::IMREAD_UNCHANGED);
        if (!rest.empty())
        {
            cv::resize(rest, rest, cv::Size(800, 400));

            cv::Mat temp[4];
            cv::split(rest, temp);
            cv::merge(temp, 4, rest);
            temp[3] *= .85;
            cv::merge(temp, 4, rest);

            rest.setTo(cv::Scalar(255, 255, 255, 255), graph);
            rest.setTo(cv::Scalar(0, 0, 0, 255), lines);
        }
    }

    cv::cvtColor(rest, rest, cv::COLOR_BGR2RGBA);
    buffer.create(rest.cols, rest.rows, rest.ptr());

    double time = 0;
    for (size_t i = 0; i < myTrajectory.size(); i++)
    {
        time += myTrajectory[i].Time;
    }
    timeNumber.setString(std::to_string(time));
    return buffer;
}


//based on https://mklimenko.github.io/english/2018/06/23/embed-resources-msvc/
class Resource
{
public:
    HGLOBAL data;
    HRSRC h;
    Resource(int id, LPWSTR type)
    {
        h = FindResource(nullptr, MAKEINTRESOURCE(id), type);
        data = LoadResource(nullptr, h);
    }

    void* loadFont(Resource fontResource)
    {
        return LockResource(fontResource.data);
    }
};


sf::Image scrollBarImageCreation()
{
    sf::Image scrollBarImage;
    cv::Mat mat(200, 15, CV_8UC3);
    mat.setTo(cv::Scalar(200, 200, 200));
    cv::line(mat, cv::Point(3, 70), cv::Point(12, 70), cv::Scalar(100, 100, 100), 1, 8);
    cv::line(mat, cv::Point(3, 100), cv::Point(12, 100), cv::Scalar(100, 100, 100), 1, 8);
    cv::line(mat, cv::Point(3, 130), cv::Point(12, 130), cv::Scalar(100, 100, 100), 1, 8);
    cv::cvtColor(mat, mat, cv::COLOR_BGR2RGBA);
    scrollBarImage.create(15, 200, mat.ptr());
    return scrollBarImage;
}

void waypointTable(std::vector<std::string> x, std::vector<std::string> y, std::vector<std::string> angle, RoundedRectangle scrollbar, RoundedRectangle textBox, sf::RenderWindow &window)
{
    std::vector<sf::VertexArray> lines;

    sf::Event event;

    sf::VertexArray xyDiv(sf::Lines, 2);
    xyDiv[0] = sf::Vertex(sf::Vector2f(58.3 + 5, 110), sf::Color(160, 160, 160));
    xyDiv[1] = sf::Vertex(sf::Vector2f(58.3 + 5, 310), sf::Color(160, 160, 160));
    lines.push_back(xyDiv);

    sf::VertexArray yaDiv(sf::Lines, 2);
    yaDiv[0] = sf::Vertex(sf::Vector2f(58.3 * 2 + 5, 110), sf::Color(160, 160, 160));
    yaDiv[1] = sf::Vertex(sf::Vector2f(58.3 * 2 + 5, 310), sf::Color(160, 160, 160));
    lines.push_back(yaDiv);

    for (size_t i = 25; i <= textBox.getLocalBounds().height; i += 25)
    {
    sf::VertexArray line25(sf::Lines, 2);
    line25[0] = sf::Vertex(sf::Vector2f(textBox.getGlobalBounds().left + 1, textBox.getGlobalBounds().top + i), sf::Color(160, 160, 160));
    line25[1] = sf::Vertex(sf::Vector2f(textBox.getGlobalBounds().left + textBox.getGlobalBounds().width - 15, textBox.getGlobalBounds().top + i), sf::Color(160, 160, 160));
    lines.push_back(line25);
    }

    for (size_t i = 0; i < lines.size(); i++)
    {
        window.draw(lines[i]);
    }
}

class Texts
{
public:
    sf::Text tabOne;
    sf::Text tabTwo;
    sf::Text x;
    sf::Text y;
    sf::Text angle;
    sf::Text time;
    sf::Text timeNumber;
    Texts(sf::Font& font)
    {
        tabOne.setFont(font);
        tabOne.setString("Graph");
        tabOne.setCharacterSize(15);
        tabOne.setFillColor(sf::Color::Black);
        tabOne.setPosition(250, 15);

        tabTwo.setFont(font);
        tabTwo.setString("Extra");
        tabTwo.setCharacterSize(15);
        tabTwo.setFillColor(sf::Color::Black);
        tabTwo.setPosition(750, 15);

        x.setFont(font);
        x.setString("x");
        x.setCharacterSize(20);
        x.setFillColor(sf::Color::Black);
        x.setPosition(50 - 20, 80);

        y.setFont(font);
        y.setString("y");
        y.setCharacterSize(20);
        y.setFillColor(sf::Color::Black);
        y.setPosition(90 - 5, 80);

        angle.setFont(font);
        angle.setString("Angle");
        angle.setCharacterSize(20);
        angle.setFillColor(sf::Color::Black);
        angle.setPosition(130, 80);

        time.setFont(font);
        time.setString("Time");
        time.setCharacterSize(20);
        time.setFillColor(sf::Color::Black);
        time.setPosition(130, 80);

        timeNumber.setFont(font);
        timeNumber.setCharacterSize(40);
        timeNumber.setFillColor(sf::Color::Black);
        timeNumber.setPosition(130, 300);
    }
};

int main()
{
    sf::RenderWindow window(sf::VideoMode(1000, 500), "AutoPilot GUI", sf::Style::Close | sf::Style::Titlebar);

    sf::RectangleShape background(sf::Vector2f((float)window.getSize().x, (float)window.getSize().y));
    background.setFillColor(sf::Color(240, 240, 240));

    sf::RectangleShape tabOne(sf::Vector2f((float)window.getSize().x / 2, 50.f));
    tabOne.setPosition(sf::Vector2f(0.f, 0.f));
    bool tabOneSelected = true;

    sf::RectangleShape tabTwo(sf::Vector2f((float)window.getSize().x / 2, 50.f));
    tabTwo.setPosition(sf::Vector2f((float)window.getSize().x / 2, 0.f));
    bool tabTwoSelected = false;

    Resource arialFont(IDR_FONT1, RT_FONT);
    sf::Font font;
    font.loadFromMemory(arialFont.loadFont(arialFont), SizeofResource(nullptr, arialFont.h));
    Texts theTexts(font);

    sf::RectangleShape graphBox(sf::Vector2f(800, 400));
    graphBox.setPosition(sf::Vector2f(200, 75));
    sf::Texture graphImage; graphImage.loadFromImage(readCVGraph("C:\\Users\\scdel\\Downloads\\57ab4a52fc09d8cf5624b2e7b34e92d5.png", theTexts.timeNumber));
    graphBox.setTexture(&graphImage);

    RoundedRectangle pointBox(5, 175, 200, 50);
    pointBox.setOutlineColor(sf::Color(180, 180, 180));
    pointBox.setOutlineThickness(1);
    pointBox.setPosition(sf::Vector2f(100, 210));

    RoundedRectangle scrollBarPlace(5, 15, 200, 50);
    scrollBarPlace.setOutlineColor(sf::Color(180, 180, 180));
    scrollBarPlace.setOutlineThickness(1);
    scrollBarPlace.setPosition(sf::Vector2f(180, 210));

    RoundedRectangle scrollBarBox(5, 15, 200, 50);
    scrollBarBox.setOutlineColor(sf::Color(180, 180, 180));
    scrollBarBox.setOutlineThickness(1);
    scrollBarBox.setPosition(sf::Vector2f(180, 210));

    sf::Texture scrollBarText;
    sf::Image scrollBarImage(scrollBarImageCreation());
    scrollBarText.loadFromImage(scrollBarImage);
    scrollBarBox.setTexture(&scrollBarText);


    std::vector<std::string> xP, yP, angles;

    sf::Cursor cursor;

    window.setFramerateLimit(30);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            switch (event.type)
            {
            case sf::Event::Closed:
                window.close();
                break;
            case sf::Event::KeyPressed:
                break;
            case sf::Event::MouseMoved:
                {
                auto x = event.mouseMove.x;
                auto y = event.mouseMove.y;
                if (tabOne.within(x, y) || tabTwo.within(x, y))
                {
                    cursor.loadFromSystem(sf::Cursor::Type::Hand);
                }
                else if (scrollBarBox.within(x, y))
                {
                    cursor.loadFromSystem(sf::Cursor::Type::Hand);
                }

                else
                {
                    cursor.loadFromSystem(sf::Cursor::Type::Arrow);
                }
                window.setMouseCursor(cursor);
                break;
                }
            case sf::Event::MouseButtonPressed:
                tabTwo.within(event.mouseButton.x, event.mouseButton.y) ? tabOneSelected = false : tabOneSelected = true;
                //if (scrollBarBox.within(event.mouseButton.x, event.mouseButton.y))
                //{
                //    double change = event.mouseButton.y;
                //    if (window.pollEvent(event) == sf::Event::MouseMoved)
                //    {
                //        change -= event.mouseButton.y;
                //        if (scrollBarBox.getPosition().y + (scrollBarBox.getGlobalBounds().height / 2) < scrollBarPlace.getGlobalBounds().height + scrollBarPlace.getGlobalBounds().top ||
                //            scrollBarBox.getPosition().y - (scrollBarBox.getGlobalBounds().height / 2) > scrollBarPlace.getGlobalBounds().height - scrollBarPlace.getGlobalBounds().top)
                //        {
                //        scrollBarBox.setPosition(scrollBarBox.getPosition().x, scrollBarBox.getPosition().y - change);
                //        }
                //        else if (scrollBarBox.getPosition().y + (scrollBarBox.getGlobalBounds().height / 2) >= scrollBarPlace.getGlobalBounds().height + scrollBarPlace.getGlobalBounds().top ||
                //            scrollBarBox.getPosition().y - (scrollBarBox.getGlobalBounds().height / 2) <= scrollBarPlace.getGlobalBounds().height - scrollBarPlace.getGlobalBounds().top)
                //        {
                //            scrollBarBox.setPosition(scrollBarBox.getGlobalBounds().height / 2 - )
                //        }
                //    }
                //}

                break;
            default:
                break;
            }
        }

        window.clear();
        window.draw(background);
        window.draw(tabOne);
        window.draw(tabTwo);
        window.draw(theTexts.tabOne);
        window.draw(theTexts.tabTwo);

        if (tabOneSelected == true)
        {
            tabOne.setFillColor(background.getFillColor());
            tabTwo.setFillColor(sf::Color(225, 225, 225));

            theTexts.tabOne.setFillColor(sf::Color::Black);
            theTexts.tabTwo.setFillColor(sf::Color(80, 80, 80));

            tabTwoSelected = false;

            window.draw(graphBox);
            window.draw(theTexts.x);
            window.draw(theTexts.y);
            window.draw(theTexts.angle);
            window.draw(theTexts.timeNumber);
            window.draw(pointBox);
            window.draw(scrollBarPlace);
            window.draw(scrollBarBox);
            //waypointTable(xP, yP, angles, scrollBarBox, pointBox, window);
        }
        else
        {
            tabTwo.setFillColor(background.getFillColor());
            tabOne.setFillColor(sf::Color(225, 225, 225));

            theTexts.tabTwo.setFillColor(sf::Color::Black);
            theTexts.tabOne.setFillColor(sf::Color(80, 80, 80));

            tabOneSelected = false;
        }

 
        window.display();
    }

    return 0;
}
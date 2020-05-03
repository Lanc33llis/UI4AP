#include "RoundedRectangle.h"


RoundedRectangle::RoundedRectangle(size_t radius, size_t width, size_t height, size_t accuracy)
{
    myRadius = radius;
    myWidth = width;
    myHeight = height;
    myAccuracy = accuracy;
    update();
}

std::size_t RoundedRectangle::getPointCount() const
{
    return myAccuracy * 4;
}

sf::Vector2f RoundedRectangle::getPoint(std::size_t index) const
{
    sf::Vector2f myPoint;
    const double pi = 3.14159265359;

    if (index <= myAccuracy)
    {
        myPoint.x = (cos(pi / 2 * index / myAccuracy) * myRadius) + (myWidth / 2 - myRadius);
        myPoint.y = (sin(pi / 2 * index / myAccuracy) * myRadius) + (myHeight / 2 - myRadius);
    }
    if (index <= myAccuracy * 2 && index > myAccuracy)
    {
        myPoint.x = (cos(pi / 2 * index / myAccuracy) * myRadius) - (myWidth / 2 - myRadius);
        myPoint.y = (sin(pi / 2 * index / myAccuracy) * myRadius) + (myHeight / 2 - myRadius);
    }
    if (index <= myAccuracy * 3 && index > myAccuracy * 2)
    {
        myPoint.x = (cos(pi / 2 * index / myAccuracy) * myRadius) - (myWidth / 2 - myRadius);
        myPoint.y = (sin(pi / 2 * index / myAccuracy) * myRadius) - (myHeight / 2 - myRadius);
    }
    if (index <= myAccuracy * 4 && index > myAccuracy * 3)
    {
        myPoint.x = (cos(pi / 2 * index / myAccuracy) * myRadius) + (myWidth / 2 - myRadius);
        myPoint.y = (sin(pi / 2 * index / myAccuracy) * myRadius) - (myHeight / 2 - myRadius);
    }
    return myPoint;
}

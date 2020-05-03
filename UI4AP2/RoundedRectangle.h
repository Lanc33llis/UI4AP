#pragma once
#include <SFML/Graphics.hpp>
class RoundedRectangle : public sf::Shape
{
    private:
        double myRadius, myWidth, myHeight;
        size_t myAccuracy;
    public:
        RoundedRectangle(size_t radius, size_t width, size_t height, size_t accuracy);
        virtual std::size_t getPointCount() const;
        virtual sf::Vector2f getPoint(std::size_t index) const;
};


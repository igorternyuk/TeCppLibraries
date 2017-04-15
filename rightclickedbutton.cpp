#include "rightclickedbutton.h"

RightClickedButton::RightClickedButton(QWidget* parent) :
    QPushButton(parent)
{}

void RightClickedButton::mousePressEvent(QMouseEvent* event)
{
    if(event->button() == Qt::LeftButton)
        emit onLeftButtonClicked();
    if(event->button() == Qt::RightButton)
        emit onRightButtonClicked();
}

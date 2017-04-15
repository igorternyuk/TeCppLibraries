#ifndef RIGHTCLICKEDBUTTON_H
#define RIGHTCLICKEDBUTTON_H

#include <QPushButton>
#include <QMouseEvent>

class RightClickedButton : public QPushButton
{
    Q_OBJECT
public:
    explicit RightClickedButton(QWidget* parent = nullptr);
signals:
    void onLeftButtonClicked();
    void onRightButtonClicked();
private slots:
    void mousePressEvent(QMouseEvent* event);

};

#endif // RIGHTCLICKEDBUTTON_H

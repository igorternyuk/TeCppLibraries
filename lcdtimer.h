#ifndef LCDTIMER_H
#define LCDTIMER_H
#include <QWidget>
#include <QLCDNumber>
#include <QTimer>
#include <QTime>

class LCDTimer : public QLCDNumber
{
    Q_OBJECT
public:
    static const int NUM_DIGIT = 8;
    static const int DISPLAY_WIDTH = 120;
    static const int DISPLAY_HEIGHT = 40;
    static const int TIMER_PERIOD = 1000;
    LCDTimer(QWidget* parent = nullptr);
    ~LCDTimer();
    const QTime* time() const;
    void reset();
private slots:
    void startTimer();
    void stopTimer();
    void update();
private:
    QTimer *timer_;
    QTime *timeElapsed_;
};

#endif // LCDTIMER_H

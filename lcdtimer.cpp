#include "lcdtimer.h"

LCDTimer::LCDTimer(QWidget *parent) :
    QLCDNumber(NUM_DIGIT, parent)
{
    timer_ = new QTimer(this);
    timeElapsed_ = new QTime(0, 0, 0);
    setFixedSize(DISPLAY_WIDTH, DISPLAY_HEIGHT);
    setAutoFillBackground(true);
    QPalette pal = this->palette();
    pal.setColor(QPalette::Normal, QPalette::WindowText, Qt::red);
    pal.setColor(QPalette::Normal, QPalette::Window, Qt::black);
    this->setPalette(pal);
    display("00:00:00");
    connect(timer_, SIGNAL(timeout()), this, SLOT(update()));
}

LCDTimer::~LCDTimer()
{
    delete timeElapsed_;
}
const QTime* LCDTimer::time() const
{
    return timeElapsed_;
}

void LCDTimer::reset()
{
    timeElapsed_->setHMS(0, 0, 0);
    display("00:00:00");
    timer_->start(TIMER_PERIOD);
}

void LCDTimer::startTimer()
{
    timer_->start(TIMER_PERIOD);
}

void LCDTimer::stopTimer()
{
    timer_->stop();
}

void LCDTimer::update()
{
    *timeElapsed_ = timeElapsed_->addMSecs(TIMER_PERIOD);
    QString timeStr = timeElapsed_->toString("hh:mm:ss");
    display(timeStr);
}

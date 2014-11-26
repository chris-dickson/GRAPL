#include <time.h>

class Time
{
public:
	Time()
	{
		Hours = 0;
		Minutes = 0;
		Seconds = 0;
	}

	Time(int seconds)
	{
		Hours	= seconds / 3600;
		Minutes = (seconds - (3600*Hours)) / 60;
		Seconds = seconds - (3600*Hours) - (60*Minutes);
	}

	Time(const Time& Other)
	{
		Hours = Other.Hours;
		Minutes = Other.Minutes;
		Seconds = Other.Seconds;
	}

	Time operator=(const Time& Other)
	{
		Hours = Other.Hours;
		Minutes = Other.Minutes;
		Seconds = Other.Seconds;

		return *this;
	}

	Time operator+(const Time& Other) const
	{
		int sec = Other.TotalSeconds() + TotalSeconds();
		return Time(sec);
	}

	Time operator-(const Time& Other) const
	{
		int sec = TotalSeconds() - Other.TotalSeconds();
		return Time(sec);
	}

	Time operator+=(const Time& Other)
	{
		Time Myself( *this );
		Myself = Myself + Other;
		Hours = Myself.Hours;
		Minutes = Myself.Minutes;
		Seconds = Myself.Seconds;
		return *this;
	}

	Time operator-=(const Time& Other)
	{
		Time Myself( *this );
		Myself = Myself - Other;
		Hours = Myself.Hours;
		Minutes = Myself.Minutes;
		Seconds = Myself.Seconds;
		return *this;
	}

	float operator/(const Time& Other) const
	{
		int MySeconds = TotalSeconds();
		int TheirSeconds = Other.TotalSeconds();

		return (float)MySeconds / (float)TheirSeconds;
	}

	int Hours;
	int Minutes;
	int Seconds;

	int TotalHours() const
	{
		return Hours;
	}

	int TotalMinutes() const
	{
		return TotalHours() * 60 + Minutes;
	}

	int TotalSeconds() const
	{
		return TotalMinutes() * 60 + Seconds;
	}
};

class GRAPLTimer
{
public:
	GRAPLTimer()
	{
		bIsTicking = false;
		bHasEverTicked = false;
	}

	void Start();
	void Stop();
	Time GetElapsed();
	void Reset();

private:

	Time		GetElapsed(time_t _start, time_t _stop);

	bool		bIsTicking;
	bool		bHasEverTicked;
	time_t		start, stop;
	Time		RunningTime;
};

void GRAPLTimer::Start()
{
	if ( bIsTicking )
	{
		return;
	}

	bHasEverTicked = true;
	bIsTicking = true;

	start = time(NULL);
}

void GRAPLTimer::Stop()
{
	if ( !bHasEverTicked )
	{
		return;
	}
	else
	{
		bIsTicking = false;
		stop = time(NULL);
		RunningTime += GetElapsed(start,stop);
	}
}


void GRAPLTimer::Reset()
{
	bIsTicking = false;
	bHasEverTicked = false;
	RunningTime.Hours = 0;
	RunningTime.Minutes = 0;
	RunningTime.Seconds = 0;

}

Time GRAPLTimer::GetElapsed()
{
	Time Ret;
	if ( bHasEverTicked )
	{
		if ( bIsTicking )
		{
			Ret = GetElapsed(start, time(NULL));
		}
		else
		{
			Ret = RunningTime;
		}
	}

	return Ret;
}

Time GRAPLTimer::GetElapsed(time_t _start, time_t _stop)
{
	Time Ret;

	double Elapsed = difftime(_stop, _start);

	Ret.Hours	= ((int)Elapsed) / 3600;
	Ret.Minutes = ((int)Elapsed - (3600*Ret.Hours)) / 60;
	Ret.Seconds = (int)Elapsed - (3600*Ret.Hours) - 60*Ret.Minutes;

	return Ret;
}


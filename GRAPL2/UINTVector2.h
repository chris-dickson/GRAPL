#pragma once
class UINTVector2
{
public:
	unsigned int x;
	unsigned int y;
};

class CompareUINTVector2
{
public:
	bool operator()(UINTVector2 s1, UINTVector2 s2)
	{
		return s1.x < s2.x || (!(s2.x < s1.x) && s1.y < s2.y);
	}
};
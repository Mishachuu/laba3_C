#pragma once
#pragma once
#pragma once
using namespace std;

class Exteption {
public:
	virtual void Print() = 0;
};


class EInvalidSize :public Exteption {
public:
	void Print();
};


class EInvalidIndex :public Exteption {
public:
	void Print();
};


class EInvalidMull :public Exteption {
public:
	void Print();
};


class EInvalidQuadrate :public Exteption {
public:
	void Print();
};


class EdivisionZero :public Exteption {
public:
	void Print();
};


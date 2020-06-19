#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <deque>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <utility>
#include <list>
#include <limits>
#include <cctype>

class Fraction
{
public:
	Fraction() : num(0), den(1) { };
	Fraction(const int64_t &n) :
		num(n), den(1) { };
	Fraction(const int64_t &n, const int64_t &d) :
		num(n), den(d) {
		this->reduece();
	};
	Fraction(const Fraction &frac) :
		num(frac.num), den(frac.den) { };
	~Fraction() { };

	friend std::ostream &operator<<(std::ostream&, const Fraction&);
	friend std::istream &operator>>(std::istream&, Fraction&); 
	friend Fraction operator+(const Fraction&, const Fraction&);
	friend Fraction operator-(const Fraction&, const Fraction&);
	friend Fraction operator*(const Fraction&, const Fraction&);
	friend Fraction operator/(const Fraction&, const Fraction&);
	friend Fraction operator/(const Fraction&, const int&);
	Fraction& operator=(const Fraction&);
	Fraction& operator=(const int&);
	Fraction operator-() const;
	friend bool operator>(const Fraction&, const int&);
	friend bool operator<(const Fraction&, const int&);
	friend bool operator>(const Fraction&, const Fraction&);
	friend bool operator<(const Fraction&, const Fraction&);
	friend bool operator>=(const Fraction&, const int&);
	friend bool operator<=(const Fraction&, const int&);
	friend bool operator>=(const Fraction&, const Fraction&);
	friend bool operator<=(const Fraction&, const Fraction&);
	friend bool operator==(const Fraction&, const int&);
	friend bool operator==(const Fraction&, const Fraction&);
	friend bool operator!=(const Fraction&, const int&);
	friend bool operator!=(const Fraction&, const Fraction&);
	Fraction& operator+=(const Fraction&);
	Fraction& operator-=(const Fraction&);
	Fraction& operator*=(const Fraction&);
	Fraction& operator/=(const Fraction&);
	Fraction& operator/=(const int&);

	inline size_t calcFracLength();
	inline const Fraction getRem() const { return Fraction(num % den, den); };

private:
	int64_t gcd(const int64_t, const int64_t);
	inline Fraction& reduece();

	int64_t num;
	int64_t den;
};

std::ostream &operator<<(std::ostream &os, const Fraction &frac)
{
	std::string fracStr = std::to_string(frac.num);
	if (frac.den != 1 && frac.num != 0)
		fracStr += '/' + std::to_string(frac.den);
	os << fracStr;
	return os;
}

std::istream &operator>>(std::istream &is, Fraction &frac)
{
	std::string inputStrFraction;
	is >> inputStrFraction;
	auto sep = inputStrFraction.find('/');

	if (sep == std::string::npos)
		frac = !std::iscntrl(inputStrFraction.front()) ? Fraction(std::stoll(inputStrFraction)) : 0;
	else
	{
		std::string strNum, strDen;
		strNum = inputStrFraction.substr(0, sep);
		strDen = inputStrFraction.substr(sep + 1, inputStrFraction.length());

		frac = Fraction(std::stoll(strNum), std::stoll(strDen));
	}

	return is;
}

Fraction& Fraction::operator+=(const Fraction &rhs)
{
	num = num * rhs.den + rhs.num * den;
	den *= rhs.den;
	this->reduece();
	return *this;
}

Fraction& Fraction::operator-=(const Fraction &rhs)
{
	num = num * rhs.den - rhs.num * den;
	den *= rhs.den;
	this->reduece();
	return *this;
}

Fraction& Fraction::operator*=(const Fraction &rhs)
{
	num *= rhs.num;
	den *= rhs.den;
	this->reduece();
	return *this;
}

Fraction& Fraction::operator/=(const Fraction &rhs)
{
	if (this == &rhs)
	{
		num = 1;
		den = 1;
		return *this;
	}
	num *= rhs.den;
	den *= rhs.num;
	this->reduece();
	return *this;
}

Fraction& Fraction::operator/=(const int &rhs)
{
	den *= rhs;
	this->reduece();
	return *this;
}

Fraction& Fraction::operator=(const Fraction &rhs)
{
	num = rhs.num;
	den = rhs.den;
	return *this;
}

Fraction& Fraction::operator=(const int &rhs)
{
	num = rhs;
	den = 1;
	return *this;
}

Fraction Fraction::operator-() const
{
	return Fraction(*this) * -1;
}

bool operator==(const Fraction &lhs, const int &rhs)
{
	return lhs.num == rhs && lhs.den == 1;
}

bool operator!=(const Fraction &lhs, const int &rhs)
{
	return !(lhs == rhs);
}

bool operator==(const Fraction &lhs, const Fraction &rhs)
{
	return lhs.num == rhs.num && lhs.den == rhs.den;
}

bool operator!=(const Fraction &lhs, const Fraction &rhs)
{
	return !(lhs == rhs);
}

bool operator>(const Fraction &lhs, const int &rhs)
{
	return lhs.num > rhs;
}

bool operator<(const Fraction &lhs, const int &rhs)
{
	return lhs.num < rhs;
}

bool operator>(const Fraction &lhs, const Fraction &rhs)
{
	return lhs.num * rhs.den > rhs.num * lhs.den;
}

bool operator<(const Fraction &lhs, const Fraction &rhs)
{
	return lhs.num * rhs.den < rhs.num * lhs.den;
}

bool operator>=(const Fraction &lhs, const int &rhs)
{
	return lhs.num >= rhs;
}

bool operator<=(const Fraction &lhs, const int &rhs)
{
	return lhs.num <= rhs;
}

bool operator>=(const Fraction &lhs, const Fraction &rhs)
{
	return lhs.num * rhs.den >= rhs.num * lhs.den;
}

bool operator<=(const Fraction &lhs, const Fraction &rhs)
{
	return lhs.num * rhs.den <= rhs.num * lhs.den;
}

Fraction operator+(const Fraction &lhs, const Fraction &rhs)
{
	Fraction sum = lhs;
	sum += rhs;
	return sum;
}

Fraction operator-(const Fraction &lhs, const Fraction &rhs)
{
	Fraction dif = lhs;
	dif -= rhs;
	return dif;
}

Fraction operator*(const Fraction &lhs, const Fraction &rhs)
{
	Fraction mult = lhs;
	mult *= rhs;
	return mult;
}

Fraction operator/(const Fraction &lhs, const Fraction &rhs)
{
	Fraction div = lhs;
	div /= rhs;
	return div;
}

Fraction operator/(const Fraction &lhs, const int &rhs)
{
	Fraction div = lhs;
	div /= rhs;
	return div;
}

inline
size_t Fraction::calcFracLength()
{
	return std::to_string(num).length() + 2 + std::to_string(den).length();
}

int64_t Fraction::gcd(const int64_t a, const int64_t b)
{
	return b == 0 ? a : gcd(b, a % b);
}

inline
Fraction& Fraction::reduece()
{
	int64_t nod = gcd(std::abs(num), std::abs(den));
	num /= nod;
	den /= nod;
	if (den < 0)
	{
		num *= -1;
		den *= -1;
	}
	return *this;
}

class SimplexTable
{
public:
	SimplexTable() : maxFracLength(0) { };
	friend std::ifstream &operator>>(std::ifstream&, SimplexTable&);
	friend std::ostream &operator<<(std::ostream&, const SimplexTable&);

	void symplexMethod();
	void methodGamori();

private:
	inline void checkMaxFracLength(Fraction&);
	inline void findBases();
	inline size_t findPivotColumnIndForSymplex();
	size_t findPivotLineIndforSymplex(const size_t&);
	size_t findPivotColumnIndForGamori();
	inline bool findNegSignInRes();
	void removeNegSignFromRes();
	void squareMethod(const size_t&, const size_t&);
	void addNewLim(const size_t&);
	inline bool checkBaseOptimality();
	inline bool checkPresenceOfNonIntegralBase();
	size_t findMinRem();
	std::deque<std::deque<Fraction>> matrix;
	std::deque<std::string> ZStr;
	std::vector<Fraction> simRel;
	size_t maxFracLength;
};

size_t SimplexTable::findPivotColumnIndForGamori()
{
	std::pair<Fraction, size_t>minelem = std::make_pair(Fraction(std::numeric_limits<int32_t>::max()), -1);
	for (size_t columnind = 0; columnind < matrix.front().size() - 2; ++columnind)
		if (matrix.front()[columnind] != 0 && minelem.first > matrix.back()[columnind] / matrix.front()[columnind])
			minelem = std::make_pair(matrix.back()[columnind] / matrix.front()[columnind], columnind);
	return minelem.second;
}

inline
bool SimplexTable::findNegSignInRes()
{
	return std::all_of(matrix.begin(), matrix.end(), [](std::deque<Fraction> &line) { return line.back() >= 0; });
}

void SimplexTable::addNewLim(const size_t &lineIndOfVar)
{
	std::deque<Fraction> newLim(matrix.front().size());
	for (size_t coulumnInd = 0; coulumnInd < matrix.front().size(); ++coulumnInd)
		newLim[coulumnInd] = matrix[lineIndOfVar][coulumnInd] >= 0 ? -matrix[lineIndOfVar][coulumnInd].getRem()
		: -matrix[lineIndOfVar][coulumnInd].getRem() - 1;
	newLim.back() = -matrix[lineIndOfVar].back().getRem();

	matrix.push_front(newLim);
	ZStr.push_back(std::string());
	std::fill(ZStr.begin(), ZStr.end(), std::string());

	for (auto &line : matrix)
		line.insert(line.end() - 1, 0);
	*(matrix.front().end() - 2) = 1;
}

void SimplexTable::removeNegSignFromRes()
{
	for (auto &line : matrix)
		if (line.back() < 0)
			for (auto &elem : line)
				elem *= -1;
}

size_t SimplexTable::findMinRem()
{
	std::vector<std::pair<Fraction, size_t>> valuseRem(matrix.size() - 1, std::make_pair(0, -1));
	for (size_t lineInd = 0; lineInd < matrix.size() - 1; ++lineInd)
	{
		size_t xInd = std::stoul(std::string().assign(ZStr[lineInd].begin() + 1, ZStr[lineInd].end())) - 1;
		if (xInd < matrix.size() - 1)
			valuseRem[xInd] = std::make_pair(matrix[lineInd].back().getRem(), lineInd);
	}
	return std::max_element(valuseRem.begin(), valuseRem.end())->second;
}

bool SimplexTable::checkPresenceOfNonIntegralBase()
{
	return std::any_of(matrix.begin(), matrix.end() - 1, [](std::deque<Fraction> &line) {return line.back().getRem() != 0; });
}

void SimplexTable::squareMethod(const size_t &baseValueY, const size_t &baseValueX)
{
	Fraction div = matrix[baseValueY][baseValueX];
	for (auto& elem : matrix[baseValueY])
		elem /= div;
	for (size_t lineInd = 0; lineInd < matrix.size(); ++lineInd)
	{
		for (size_t columnInd = 0; columnInd < matrix[lineInd].size(); ++columnInd)
			if (lineInd != baseValueY && columnInd != baseValueX)
			{
				matrix[lineInd][columnInd] = matrix[baseValueY][baseValueX] * matrix[lineInd][columnInd] -
					matrix[baseValueY][columnInd] * matrix[lineInd][baseValueX];

				checkMaxFracLength(matrix[lineInd][columnInd]);
			}

		if (lineInd != baseValueY)
			matrix[lineInd][baseValueX] = 0;
	}
}

inline
bool SimplexTable::checkBaseOptimality()
{
	return std::all_of(matrix.back().begin(), matrix.back().end(), [](const Fraction &it) { return it >= 0; });
}

inline
size_t SimplexTable::findPivotColumnIndForSymplex()
{
	return std::min_element(std::begin(matrix.back()), std::end(matrix.back())) - matrix.back().begin();
}

size_t SimplexTable::findPivotLineIndforSymplex(const size_t &pivotColumnInd)
{
	std::pair<Fraction, size_t> minElem = std::make_pair(Fraction(std::numeric_limits<int32_t>::max()), -1);
	for (size_t lineInd = 0; lineInd < matrix.size() - 1; ++lineInd)
	{
		Fraction frac = matrix[lineInd][pivotColumnInd] > 0 ? matrix[lineInd].back() / matrix[lineInd][pivotColumnInd] : -1;
		simRel.push_back(frac);
		if (frac >= 0 && frac < minElem.first)
			minElem = std::make_pair(frac, lineInd);
	}
	return minElem.second;
}

inline
void SimplexTable::checkMaxFracLength(Fraction &frac)
{
	size_t fracLength = frac.calcFracLength();
	maxFracLength = fracLength > maxFracLength ?
		fracLength : maxFracLength;
}

std::ostream &operator<<(std::ostream &os, const SimplexTable &table)
{
	os << std::setw(table.maxFracLength + 1) << "bases";
	for (size_t xInd = 0; xInd < table.matrix.front().size() - 1; ++xInd)
		os << std::setw(table.maxFracLength + 1) << "x" + std::to_string(xInd + 1);
	os << std::setw(table.maxFracLength + 1) << "eq" << std::setw(table.maxFracLength + 1) << "SR" << '\n';
	for (size_t lineInd = 0; lineInd < table.matrix.size(); ++lineInd)
	{
		os << std::setw(table.maxFracLength + 1) << table.ZStr[lineInd];
		for (auto frac : table.matrix[lineInd])
			os << std::setw(table.maxFracLength + 1) << frac;
		os << std::setw(table.maxFracLength + 1);
		try
		{
			if (table.simRel.at(lineInd) >= 0)
				os << table.simRel.at(lineInd);
			else
				os << '-';
		}
		catch (std::out_of_range)
		{
			os << "N\\C";
		}
		os << '\n';
	}
	return os;
}


void SimplexTable::findBases()
{
	std::fill(ZStr.begin(), ZStr.end(), "");

	for (size_t columnInd = 0; columnInd < matrix.front().size(); ++columnInd)
	{
		uint16_t zeroesAmount = 0;
		size_t baseVarInd;

		for (size_t lineInd = 0; lineInd < matrix.size(); ++lineInd)
			matrix[lineInd][columnInd] == 0 ? ++zeroesAmount : baseVarInd = lineInd;

		if (zeroesAmount == matrix.size() - 1 && matrix[baseVarInd][columnInd] == 1)
			ZStr[baseVarInd] = 'x' + std::to_string(columnInd + 1);
	}

	for (size_t lineInd = 0; lineInd < ZStr.size() - 1; ++lineInd)
		if (ZStr[lineInd] == "")
		{
			size_t columnInd = std::find_if(matrix[lineInd].begin(), matrix[lineInd].end(),
				[](const Fraction &it) {return it > 0; }) - matrix[lineInd].begin();
			ZStr[lineInd] = 'x' + std::to_string(columnInd + 1);
			squareMethod(lineInd, columnInd);
		}

	ZStr.back() = 'Z';
}

std::ifstream &operator>>(std::ifstream &is, SimplexTable &table)
{
	table.matrix.clear();
	std::string EqFromFile;

	while (std::getline(is, EqFromFile))
	{
		std::stringstream strEq(EqFromFile);

		Fraction frac;
		std::deque<Fraction> Eq;
		while (strEq >> frac)
		{
			Eq.push_back(frac);
			table.checkMaxFracLength(frac);
		}

		table.matrix.push_front(Eq);
	}

	for (size_t i = 0; i < table.matrix.size() - 1; ++i)
	{
		Fraction eqRes = table.matrix[i].back();
		table.matrix[i].pop_back();
		for (size_t j = 0; j < table.matrix.size() - 1; ++j)
			table.matrix[i].push_back(i == j ? 1 : 0);
		table.matrix[i].push_back(eqRes);
	}

	for (auto &elem : table.matrix.back())
		elem *= -1;

	while (table.matrix.back().size() != table.matrix.front().size())
		table.matrix.back().push_back(0);

	table.ZStr.resize(table.matrix.size());
	table.findBases();
}

void SimplexTable::symplexMethod()
{
	while (!checkBaseOptimality())
	{
		std::cout << *this << '\n';
		size_t pivotColumnInd = findPivotColumnIndForSymplex();
		size_t pivotLineInd = findPivotLineIndforSymplex(pivotColumnInd);
		std::cout << *this << '\n';
		while (!findNegSignInRes())
		{
			removeNegSignFromRes();
			squareMethod(pivotLineInd, pivotColumnInd);
		}
		squareMethod(pivotLineInd, pivotColumnInd);
		simRel.clear();
		findBases();
	}
	std::cout << *this;
}

void SimplexTable::methodGamori()
{
	this->symplexMethod();
	if (checkPresenceOfNonIntegralBase())
		std::for_each(matrix.back().begin(), matrix.back().end(), [](Fraction &res) { res = -res; });
	while (checkPresenceOfNonIntegralBase())
	{
		addNewLim(findMinRem());
		findBases();
		std::cout << '\n' << *this;
		size_t pivotColumn = findPivotColumnIndForGamori();
		squareMethod(0, pivotColumn);
		findBases();
		std::cout << '\n' << *this;
	}
}

int main()
{
	SimplexTable table;
	std::ifstream in("test5.txt");
	in >> table;
	table.methodGamori();
}

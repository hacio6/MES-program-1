#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

struct wezel {

	double x;		//wspó³rzêdna
	int WB;		//warunki brzegowe
	wezel(int x, int WB) {
		this->x = x;
		this->WB = WB;	
	}
	wezel() {
		this->x = 0;
		this->WB = 0;
	}
};
struct element {

	int nr_el;	//numer elementu
	int id1, id2;	//wêze³ pocz¹tkowy i wêze³ koñcowy
	double K;	//wspó³czynnik przewodzenia ciep³a
	double S;	//powierzchnia na której zadane s¹ warunki brzegowe
	double L;	//d³ugoœæ 
	double **lH;	//lokalna macierz H
	double *lP;		// lokalny wektor P
	wezel *w;
	element(int nr_el, int id1, int id2) {
		this->nr_el = nr_el;
		this->id1 = id1;
		this->id2 = id2;
		this->w = new wezel[2];
	}
	element() {
		this->id1 = 0;
		this->id2 = 0;
		this->nr_el = 0;
		this->w = new wezel[2];
	}
};
struct uklad {
	double **gH;	//globalna macierz H
	double *gP;		//globalny wektor P
	bool gauss(int n, double **AB, double *X) // metoda algorytmu eliminacji Gaussa
	{
		const double eps = 1e-12; // sta³a przybli¿enia zera
		int i, j, k;
		double m, s;

		// eliminacja wspó³czynników

		for (i = 0; i < n - 1; i++)
		{
			for (j = i + 1; j < n; j++)
			{
				if (fabs(AB[i][i]) < eps) return false;
				m = -AB[j][i] / AB[i][i];
				for (k = i + 1; k <= n; k++)
					AB[j][k] += m * AB[i][k];
			}
		}

		// wyliczanie niewiadomych

		for (i = n - 1; i >= 0; i--)
		{
			s = AB[i][n];
			for (j = n - 1; j >= i + 1; j--)
				s -= AB[i][j] * X[j];
			if (fabs(AB[i][i]) < eps) return false;
			X[i] = s / AB[i][i];
		}
		return true;
	}
};
struct siatka {
	int ne;		//liczba elementów
	int nh;		//liczba wez³ów
	double q, alfa, tot;	//strumieñ ciep³a q, wspó³czynnik wymiany ciep³a, temperatura otoczenia
	element *el;
	uklad ur;
	siatka(int ne) {
		this->ne = ne;
		this->nh = ne + 1;
		el = new element[ne];
	}
	siatka() {
		this->ne = 0;
		this->nh = 0;
		el = new element[ne];
	}
};


int main() {
	/*-------------------------------------------------------------------------------------
	---------------------------WCZYTYWANIE DANYCH Z PLIKU--------------------------------
	-------------------------------------------------------------------------------------*/
	int n;
	double x = 0;	//zmienna pomocnicza do okreœlenia wspó³rzêdnej 
	fstream p, wynik;
	p.open("daneW.txt", ios::in | ios::out);
	if (p.good() == true) {
		p >> n;
		siatka *siat = new siatka(n);
		for (int i = 0; i < n; i++)
		{

			int j = 0;
			p >> siat->el[i].id1;	//id wezla pocz
			p >> siat->el[i].w[j].WB;	//warunek brzegowy pocz		
			p >> siat->el[i].id2;	//id wezla kon
			siat->el[i].w[j].x = x;
			j++;
			p >> siat->el[i].w[j].WB;	//warunek brzegowy wspolnego wezla
			p >> siat->el[i].K; //wspolczynnik przewodzenia ciepla w tym elemencie
			p >> siat->el[i].S;	//powierzchnia elementu
			p >> siat->el[i].L;	//dlugosc elementu
			siat->el[i].w[j].x = x + siat->el[i].L;

			cout << "Dane elementu " << i + 1 << endl;
			cout << "Wezel poczatkowy = " << siat->el[i].id1 << "	Wezel koncowy = " << siat->el[i].id2 << "	K = " << siat->el[i].K << "	S = " << siat->el[i].S << " L = " << siat->el[i].L << endl;
			cout << "Warunek brzegowy = " << siat->el[i].w[j - 1].WB << "	Warunek brzegowy z drugiej strony = " << siat->el[i].w[j].WB << "	Wspolrzedna wezla pocz = " << siat->el[i].w[j-1].x << "	Wspolrzedna wezla kon = " << siat->el[i].w[j].x << endl;
			cout << endl;

		}
		p >> siat->q;
		p >> siat->alfa;
		p >> siat->tot;
		cout << "q = " << siat->q << "	alfa = " << siat->alfa << "	tot = " << siat->tot << endl;
		cout << endl;
		/*-------------------------------------------------------------------------------------
		--------------------------------------MACIERZE---------------------------------------
		-------------------------------------------------------------------------------------*/
		//		MACIERZ ALFA*S
		double macierzAlfaS[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				macierzAlfaS[i][j] = 0;
		macierzAlfaS[1][1] = siat->alfa * siat->el[0].S;
		//		MACIERZE LOKALNE H
		cout << "MACIERZE LOKALNE H" << endl;
		cout << endl;
		for (int i = 0; i < n; i++)
		{
			siat->el[i].lH = new double *[n];
			for (int j = 0; j < 2; j++)
			{
				siat->el[i].lH[j] = new double[n];
			}

		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					if (j == k) {
						siat->el[i].lH[j][k] = (siat->el[i].K*siat->el[i].S) / siat->el[i].L;
					}
					else
					{
						siat->el[i].lH[j][k] = -(siat->el[i].K*siat->el[i].S) / siat->el[i].L;
					}
					if (siat->el[i].w[j].WB == 2)
					{
						siat->el[i].lH[j][k] += macierzAlfaS[j][k];
					}
				}
			}

		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					cout << siat->el[i].lH[j][k] << "	";
				}
				cout << endl;
			}
			cout << endl;
		}

		//GLOBALNA MACIERZ H
		cout << "MACIERZ GLOBALNA H" << endl;
		cout << endl;

		siat->ur.gH = new double *[n + 1];
		for (int i = 0; i < (n+1); i++)
		{
			siat->ur.gH[i] = new double[n + 1];
		}

		for (int i = 0; i < (n + 1); i++)
			for (int j = 0; j < (n + 1); j++)
				siat->ur.gH[i][j] = 0;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					siat->ur.gH[j + i][k + i] += siat->el[i].lH[j][k];
				}
			}
		}

		for (int j = 0; j < (n+1); j++)
		{
			for (int k = 0; k < (n+1); k++)
			{
				cout << siat->ur.gH[j][k] << "	";
			}
			cout << endl;
		}
		cout << endl;

		//WEKTORY LOKALNE P
		cout << "WEKTORY LOKALNE P" << endl;
		cout << endl;
		for (int i = 0; i < n; i++)
		{
			siat->el[i].lP = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (siat->el[i].w[j].WB == 0)
					siat->el[i].lP[j] = 0;
				else if (siat->el[i].w[j].WB == 1)
					siat->el[i].lP[j] = siat->q*siat->el[i].S;
				else if (siat->el[i].w[j].WB == 2)
					siat->el[i].lP[j] = -(siat->alfa*siat->tot*siat->el[i].S);
			}
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				cout << siat->el[i].lP[j] << "	" << endl;
			}
			cout << endl;
		}
		cout << endl;

		//WEKTOR GLOBALNY P
		cout << "WEKTOR GLOBALNY P" << endl;
		cout << endl;
		for (int i = 0; i < n; i++)
		{
			siat->ur.gP = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 1; k < n; k++) {
					if (siat->el[i].w[j].WB == 0)
						siat->ur.gP[k] = 0;
					else if (siat->el[i].w[j].WB == 1)
						siat->ur.gP[0] = siat->q*siat->el[i].S;
					else if (siat->el[i].w[j].WB == 2)
						siat->ur.gP[n] = -(siat->alfa*siat->tot*siat->el[i].S);
				}
			}
		}

		for (int i = 0; i < (n + 1); i++)
		{
			cout << siat->ur.gP[i] << "	" << endl;
		}
		cout << endl;




		double **HP;
		HP = new double *[n + 1];
		for (int i = 0; i < (n + 1); i++)
			HP[i] = new double[n + 2];

		for (int i = 0; i < (n + 1); i++)
			for (int j = 0; j <= (n + 1); j++)
				HP[i][j] = 0;


		for (int i = 0; i < (n + 1); i++)
		{
			for (int j = 0; j <= (n + 1); j++)
			{
				HP[i][j] = siat->ur.gH[i][j];
				HP[i][n + 1] = -siat->ur.gP[i];
			}
		}



		double *t;
		t = new double[n + 1];
		cout << "WYNIKI" << endl;
		cout << endl;
		wynik.open("wynikiW.txt", ios::in | ios::out);
		if (wynik.good() == true) {
			if (siat->ur.gauss((n + 1), HP, t))
			{
				for (int i = 0; i < (n + 1); i++) {
					cout << "t" << i + 1 << " = " << t[i] << endl;
					wynik << t[i] << endl;
				}
			}
			else
				cout << "COS NIE TAK!\n";
			wynik.close();
		}
		else cout << "Brak dostepu do pliku!" << endl;
		p.close();
	}
	else cout << "Brak dostepu do pliku!" << endl;




	system("pause");
	return 0;
}
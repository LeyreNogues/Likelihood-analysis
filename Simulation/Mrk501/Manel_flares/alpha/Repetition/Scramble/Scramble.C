void Scramble(char *fname){

	const Int_t numberOfEvents = 1491;
	Double_t E[numberOfEvents],t[numberOfEvents];

	ifstream in;
	in.open(fname);
	double v1,v2;
	int Npoints = 0;
	while(1)
	{
	  in >> v1 >> v2;
	  if (!in.good()) break;
	  t[Npoints] = v1;
	  E[Npoints] = v2;
	  Npoints++;
	}

	//Scramble times
	for(int i=numberOfEvents-1; i>0; i--){
			int ord = rand()%i;
			Double_t temp = t[i];
			t[i]=t[ord];
			t[ord]=temp;
		}

	//Create a txt file with the events
		FILE *fout;
		fout = fopen("scramble_test.txt","wb");
		for(int i = 0; i<numberOfEvents; i++)
		{
			fprintf(fout,"%0.8lf    %0.8lf   \n", t[i], E[i]);
		}

}

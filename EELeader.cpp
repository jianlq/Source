#include"evolutionbit.h"
#include"EE.h"

//stackelberg4 能耗和吞吐率
int main(){
	srand((unsigned)time(NULL));
	int Time = 3;
	int CASEnum= 30;	
	vector<double> STARTUP;
	STARTUP.push_back(0);
	STARTUP.push_back(10);
	STARTUP.push_back(100);
	STARTUP.push_back(500);
	STARTUP.push_back(1000);
	STARTUP.push_back(2000);
	STARTUP.push_back(4000);
	STARTUP.push_back(10000);

	vector<double>CONSIDER;
	int CON_VALUE = 25;
	double c =0;
	for(int i=0;i <= CON_VALUE;i++){
		CONSIDER.push_back(c);
		c += 0.2;
	}

	for(int i =0;i<Time;i++){
		for(int start = 0; start < STARTUP.size();start++){
			FILE *out = fopen("outputFile//eeor.csv", "a");
			FILE *res = fopen("outputFile//result.csv", "a");


			int conN = CONSIDER.size();
			vector<double> cmpmlu(conN,0) ;
			vector<double> cmpdelay(conN,0);
			vector<int> successCase (conN, 0) ;
			vector<int> flag(conN,1);

			for(int casenum = 0; casenum < CASEnum; casenum++){
				genGraph(9,50,"inputFile//graph.txt");
				genGraphOR(9,4,9,"inputFile//graphOR.txt");
				CGraph *G = new CGraph("inputFile//graph.txt");
				CGraph *GOR = new CGraph("inputFile//graphOR.txt");

				vector<demand> eqOR;
				vector<demand> eqTE;
				vector<demand> eqbase;

				//eqbase.clear();//background流 
				for(int i = 0; i < BGNUM; i++){
					int s = rand()%G->n, t;
					do{
						t = rand()%G->n;
					}while( s == t || G->canNotReach(s,t));
					eqbase.push_back(demand(s, t, rand()%(MAXDEMAND)+1));
				}

				////Overlay  产生demand
				//eqOR.clear(); 
				for(int i =0 ;i<GOR->m;i++)
					if(G->canNotReach(GOR->Link[i]->tail, GOR->Link[i]->head))
						continue;
					else
						eqOR.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,rand()%(MAXDEMAND)+1));

				//eqTE.clear(); 
				for(unsigned int i=0;i<eqOR.size();i++)
					eqTE.push_back(eqOR[i]);
				int ornum = eqTE.size();
				for(unsigned int i=0;i<eqbase.size();i++)
					eqTE.push_back(eqbase[i]);

				double eedic = 0, ordic = 0;

				G->clearOcc();
				eedic = EEdictor(G,eqTE,ornum,STARTUP[start]);

				G->clearOcc();
				ordic = throughput(G,eqTE,ornum,STARTUP[start]);

				exit(1);

				G->clearOcc();
				if(!G->GAinit(eqTE)){
					cout << "*****GA init failed***** "<<endl;
					break;
				}

				if( (eedic + 1e-5 >= INF)  || ( ordic + 1e-5 >=INF) )
					break;

				for(int con = 0;con < CONSIDER.size();con++){


					int n = 100;//种群个体数目
					int m = eqTE.size();
					evoluPopubit popubit(n,m,G,GOR,&eqTE,&eqOR,eedic,ordic,CONSIDER[con],STARTUP[start]);
					(popubit).evolution();
					cout<<"S\t"<<popubit.hero.energy<<"\t"<<popubit.hero.delay<<endl;

					if( (popubit.hero.energy +1e-5) >= INF ||  (popubit.hero.delay + 1e-5) > INF ){
						flag[con] = 0;
						break;
					}

					if(flag[con]){
						successCase[con] += 1;
						cmpmlu[con] += popubit.hero.energy/eedic;
						cmpdelay[con] += popubit.hero.delay/ordic;	
						fprintf(out,"EE,%f,%f,%f\n",STARTUP[start],eedic,G->delay);
						fprintf(out,"OR,%f,%f,%f\n",STARTUP[start],G->energy,ordic);
						fprintf(out,"S,%f,%f,%f\n",STARTUP[start],popubit.hero.energy,popubit.hero.delay);
						fclose(out);

					}
				}
				delete G;
				delete GOR;
			}
			for(int con = 0;con < CONSIDER.size();con++)
				fprintf(res, "%f,%f,%f,%f\n",STARTUP[start],CONSIDER[con],cmpmlu[con]/successCase[con],cmpdelay[con]/successCase[con]); 		
			fclose(res);
		}
	}
	system("pause");
	return 0;	
}
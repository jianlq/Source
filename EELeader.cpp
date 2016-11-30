#include"evolutionbit.h"
#include"EE.h"
#include"nash.h"

// 能耗和吞吐率
int main(){
	srand((unsigned)time(NULL));
	int Time = 3;
	int CASEnum= 20;	
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
	double c = 0;
	for(int i=0;i <= CON_VALUE;i++){
		CONSIDER.push_back(c);
		c += 0.2;
	}
	int LOOP  = 100;

	for(int i =0;i<Time;i++){
		for(unsigned int start = 0; start < STARTUP.size();start++){
			int conN = CONSIDER.size();
			vector<double> see(conN,0) ;
			vector<double> sbw(conN,0);

			vector<double> ee(conN,0) ;
			vector<double> eebw(conN,0) ;

			vector<double> bwee(conN,0);
			vector<double> bw(conN,0);

			vector<double> nashee(conN,0) ;
			vector<double> nashbw(conN,0);

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

				G->clearOcc();
				if(!G->GAinit(eqTE)){
					cout << "*****GA init failed***** "<<endl;
					break;
				}

				if( (eedic + 1e-5 >= INF)  || ( ordic - 1e-5 <=SMALL) )
					break;

				//// nash	
				FILE *nash = fopen("outputFile//nash.csv", "a");
				int nacase = 0;
				double loopnashee=0,loopnashbw=0;
				fprintf(nash,"\n\n nash \n");
				fprintf(nash,"STARTUP,%f \n",STARTUP[start]);
				for(int i =0;i<LOOP;i++){
					G->clearOcc();
					GOR->clearOcc();
					double curee = NashEE(G,GOR,eqTE,STARTUP[start]);
					if(curee + 1e-5 >= INF){
						fprintf(nash,"NashEE unfeasible\n");
						break;
					}
					double curbw = bwcplex(GOR,eqOR);
					if( curbw - 1e-5 <= SMALL){
						fprintf(nash,"NashOR unfeasible\n");
						break;
					}
					eqTE.clear();
					for(int i=0;i<GOR->m;i++){
						if(GOR->Link[i]->use>0)
							eqTE.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,GOR->Link[i]->use));
					}
					for(unsigned int i=0;i<eqbase.size();i++)
						eqTE.push_back(eqbase[i]);

					nacase++;
					loopnashee += curee;
					loopnashbw += curbw;
					fprintf(nash,"%f,%f\n",curee,curbw);
					printf("%f,%f\n",curee,curbw);
				}
				fclose(nash);

				FILE *cur = fopen("outputFile//eebw.csv", "a");
				fprintf(cur,"%f\n",STARTUP[start]);
				fprintf(cur,"\n\n,EE,,%f,%f\n",eedic,G->throughput);
				fprintf(cur,",OR,,%f,%f\n",G->energy,ordic);
				fprintf(cur,",Nash,,%f,%f\n\n",loopnashee/nacase,loopnashbw/nacase);
				fclose(cur);

				for(unsigned int con = 0;con < CONSIDER.size();con++){

					int n = 150;//种群个体数目
					int m = eqTE.size();
					evoluPopubit popubit(n,m,G,GOR,&eqTE,&eqOR,eedic,ordic,CONSIDER[con],STARTUP[start]);
					(popubit).evolution();
					cout<<"S\t"<<popubit.hero.energy<<"\t"<<popubit.hero.throughput <<endl;

					if( (popubit.hero.energy +1e-5) >= INF ||  (popubit.hero.throughput  - 1e-5) < SMALL ){
						flag[con] = 0;
						break;
					}

					cur = fopen("outputFile//eebw.csv", "a");
					if(flag[con]){	
						fprintf(cur,",S,%f,%f,%f\n",CONSIDER[con],popubit.hero.energy,popubit.hero.throughput);
						fclose(cur);

						successCase[con] += 1;

						ee[con] += eedic;
						eebw[con] += G->throughput;

						bwee[con] += G->energy;
						bw[con] += ordic;						

						see[con] += popubit.hero.energy;
						sbw[con] += popubit.hero.throughput;

						nashee[con] += (loopnashee/nacase);
						nashbw[con] += (loopnashbw/nacase);
						
					}

				} // end of CONSIDER for

				delete G;
				delete GOR;

			} // end of CASENum for

			FILE *res = fopen("outputFile//result.csv", "a");		
			fprintf(res, "%f\n",STARTUP[start]);
			fprintf(res,"%d case average\n",CASEnum);
			fprintf(res,"successCase,CONSIDER,Energy Efficiency,,,,throughput,,,,Stackelberg,,,,Nash\n");
			for(unsigned int con = 0;con < CONSIDER.size();con++){
				fprintf(res, "EE,%f,%f,%f,,OR,%f,%f,,S,%f,%f,,nash,%f,%f\n",CONSIDER[con],ee[con]/successCase[con],eebw[con]/successCase[con]
				,bw[con]/successCase[con],bwee[con]/successCase[con]
				,see[con]/successCase[con],sbw[con]/successCase[con]
				,nashee[con]/successCase[con],nashbw[con]/successCase[con]); 
			}
			fclose(res);

		} // end of STARTUP for

	} // end of Time for
	system("pause");
	return 0;	
}
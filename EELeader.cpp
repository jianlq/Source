#include"evolutionbit.h"
#include"EE.h"
#include"nash.h"

// 能耗和吞吐率
int main(){
	srand((unsigned)time(NULL));
	int Time = 2;
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

			double ee = 0;
			double eebw = 0;

			double bwee = 0;
			double bw = 0;

			double nashee = 0;
			double nashbw = 0;

			vector<int> successCase (conN, 0) ;

			int sucCaseEE = 0,sucCaseBW = 0,sucCaseNash = 0;

			for(int casenum = 0; casenum < CASEnum; casenum++){

				genGraph(9,54,"inputFile//graph.txt");
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

				if( eedic + 1e-5 <= INF ){
					sucCaseEE++;
					ee += eedic;
					eebw += G->throughput;
				}
				else
					break;

				G->clearOcc();
				ordic = throughput(G,eqTE,ornum,STARTUP[start]);

				if( ordic - 1e-5 >=SMALL ){
					sucCaseBW++;
					bw += ordic;
					bwee += G->energy;
				}
				else
					break;
					
				G->clearOcc();
				if(!G->GAinit(eqTE)){
					cout << "*****GA init failed***** "<<endl;
					break;
				}

				// 让两个dic打？nash 结果？？？ 2016-12-01 12:00

				//// nash		
				FILE *nash = fopen("outputFile//nash.csv", "a");
				int nacase = 0;
				double loopnashee=0,loopnashbw=0;
				fprintf(nash,"\n\n nash \n");
				fprintf(nash,"STARTUP,%d \n",STARTUP[start]);
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
					fprintf(nash,"%lf,%lf\n",curee,curbw);
					cout<<curee<<"\t"<<curbw<<endl;
				}
				fclose(nash);

				if(nacase){
					nashee += (loopnashee/nacase);
					nashbw += (loopnashbw/nacase);
					sucCaseNash++;
				}

				FILE *cur = fopen("outputFile//eebw.csv", "a");
				fprintf(cur,"\n\n\n%lf,%d\n",STARTUP[start],casenum);
				fprintf(cur,",EE,,%lf,%lf\n",eedic,G->throughput);
				fprintf(cur,",OR,,%lf,%lf\n",G->energy,ordic);
				fprintf(cur,",Nash,,%lf,%lf\n\n",loopnashee/nacase,loopnashbw/nacase);
				fclose(cur);

				for(unsigned int con = 0;con < CONSIDER.size();con++){

					int n = 150;//种群个体数目
					int m = eqTE.size();
					evoluPopubit popubit(n,m,G,GOR,&eqTE,&eqOR,eedic,ordic,CONSIDER[con],STARTUP[start]);
					(popubit).evolution();
					cout<<"S\t"<<popubit.hero.energy<<"\t"<<popubit.hero.throughput <<endl;

					if( (popubit.hero.energy +1e-5) >= INF ||  (popubit.hero.throughput  - 1e-5) < SMALL ){
						break;
					}

					cur = fopen("outputFile//eebw.csv", "a");
						
					fprintf(cur,",S,%lf,%lf,%lf\n",CONSIDER[con],popubit.hero.energy,popubit.hero.throughput);
					fclose(cur);

					successCase[con] += 1;
					see[con] += popubit.hero.energy;
					sbw[con] += popubit.hero.throughput;

				} // end of CONSIDER for

				delete G;
				delete GOR;

			} // end of CASENum for

			FILE *res = fopen("outputFile//result.csv", "a");		
			fprintf(res, "\n\nSTARTUP,%lf\n",STARTUP[start]);
			fprintf(res,"%d case average\n",CASEnum);
			fprintf(res,",,CONSIDER,successCase,Energy Efficiency,throughput\n");
			for(unsigned int con = 0;con < CONSIDER.size();con++){
				fprintf(res, ",S,%lf,%d,%lf,%lf\n",CONSIDER[con],successCase[con],see[con]/successCase[con],sbw[con]/successCase[con]); 
			}
			fprintf(res,",,,successCase,Energy Efficiency,throughput\n");
			fprintf(res, "\n\n,EE,,%d,%lf,%lf\n",sucCaseEE,ee/sucCaseEE,eebw/sucCaseEE);
			fprintf(res, ",OR,,%d,%lf,%lf\n",sucCaseBW,bw/sucCaseBW,bwee/sucCaseBW);
			fprintf(res, ",Nash,,%d,%lf,%lf\n",sucCaseNash,nashee/sucCaseNash,nashbw/sucCaseNash);
			fclose(res);

		} // end of STARTUP for

	} // end of Time for
	system("pause");
	return 0;	
}
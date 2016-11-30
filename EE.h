#ifndef EE_H
#define EE_H
#include"CGraph.h"
#include <ilcplex/ilocplex.h>

// OPEN+x^2
double EEdictor(CGraph *G,vector<demand> & req,int ornum,double OPEN){
	IloEnv env;
	IloModel model(env);
	IloCplex EEsolver(model);

	int num = req.size();	
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); 	

	// 对每个点进行流量守恒约束  
	for(int d = 0; d < num; d++){
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];

			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	for(int i = 0; i < G->m; i++){
		IloExpr constraint(env);
		for(int d = 0; d <  num; d++)
			constraint += req[d].flow*x[d][i];
		model.add( constraint <= G->Link[i]->capacity );  
	}

	//优化目标 节能 OPEN+X^2
	IloExpr cost(env);
	for(int i = 0; i < G->m; i++){
		IloExpr load(env);
		IloIntVarArray tmp(env,num,0,1);
		for(int d = 0; d < num; d++){
			load += req[d].flow*x[d][i];
			tmp[d] = x[d][i];
		}
		//cost += ( IloPower(load,2) + IloMax(tmp)*OPEN );
		cost += ( load*load + IloMax(tmp)*OPEN );
	}
	model.add(IloMinimize(env,cost));

	EEsolver.setOut(env.getNullStream());
	double obj = INF;
	if(EEsolver.solve()){
		obj = EEsolver.getObjValue(); //energy
		//thoughtput
		double output = 0;
		for(int i = 0;i < G->m; i++){  
			double loadc = 0;
			for(int d = 0;d < num; d++)
				loadc += EEsolver.getValue(x[d][i])*req[d].flow;
			G->Link[i]->bw = G->Link[i]->capacity -loadc; //residual bw
		}
		for(int d = 0;d < ornum; d++){
			double res = INF;
			for(int i = 0;i < G->m; i++){
				if(EEsolver.getValue(x[d][i]) > 0.5)
					res = min( res, G->Link[i]->bw );
			}
			output += res*req[d].flow;
		}
		G->throughput = output;
		cout << "EE\t能耗 "<< obj <<"\t吞吐率 "<<output<<endl;
	}
	else{
		cout << "EE unfeasible"<<endl;
	}
	for(size_t i = 0; i < req.size(); i++)
		x[i].end();
	x.end();
	env.end();
	return obj;
}

//数据流 吞吐率 P18
double throughput(CGraph *G,vector<demand> req,int ornum,double OPEN){    
	IloEnv env;
	IloModel model(env);
	IloCplex ORsolver(model);

	int totalnum = req.size();
	IloArray<IloIntVarArray> x(env, totalnum); 
	for(int d = 0; d < totalnum; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); // num * G->m 的二维矩阵

	//优化目标
	IloExpr goal(env);
	IloExprArray Load(env, G->m);
	IloIntVarArray Res(env,totalnum,0,INF); //Res[d]表示业务d可用的剩余带宽
	// IloInfinity无解 ？？？   INF

	for(int i=0;i<G->m;i++)
		Load[i] = IloExpr(env);


	for(int i = 0; i < G->m; i++){
		for(int d = 0; d < totalnum; d++)
			Load[i] += x[d][i]*req[d].flow;
	}

	for(int d = 0; d < ornum; d++){
		for(int i = 0; i < G->m; i++)
			model.add(Res[d] <= ((1-x[d][i])*INF + (G->Link[i]->capacity - Load[i])));
		goal += Res[d]*req[d].flow; //why req[d].flow
		//goal+=Y[d]*reqL[d]->flow/(int)g->cost_best[reqL[d]->id];
	}

	model.add(IloMaximize(env,goal));

	//流量守恒约束
	for(int d = 0; d < totalnum; d++)
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];
			// 出 - 入
			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}

		//带宽约束
		for(int i = 0; i < G->m; i++){
			IloExpr constraint(env);
			for(int d = 0; d <  totalnum; d++)
				constraint += req[d].flow*x[d][i];
			model.add(constraint <= G->Link[i]->capacity);  
		}

		ORsolver.setOut(env.getNullStream());
		double obj = SMALL;
		ORsolver.solve();
		if(ORsolver.getStatus() == IloAlgorithm::Infeasible)
			env.out() << "throughput unfeasible" << endl;
		else{
			obj = ORsolver.getObjValue();
			// 计算 energy
			double energy = 0;
			for(int i = 0; i < G->m; i++){  
				double load = 0;
				double one = 0;
				for(int d = 0; d < totalnum; d++){
					load += ORsolver.getValue(x[d][i])*req[d].flow;
					one = max(one,ORsolver.getValue(x[d][i]));
				}
				energy += (one*OPEN + load*load);
			}
			G->energy = energy;
			cout << "OR\t能耗 "<<energy<< "\t吞吐率 "<<obj<<endl;
		}
		for(size_t i = 0; i < req.size(); i++)
			x[i].end();
		x.end();
		env.end();
		return obj;
}

//nashbw
double bwcplex(CGraph *g,vector<demand> req){    
	IloEnv env;
	IloModel model(env);
	IloCplex ORsolver(model);

	int num = req.size();
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, g->m, 0, 1); 

	//优化目标
	IloExpr goal(env);
	
	IloIntVarArray cost(env,num,0,INF); //cost[d]表示业务d可用的剩余带宽

	for(int d = 0; d < num; d++){
		for(int i = 0; i < g->m; i++)
			model.add(cost[d] <= ( (1-x[d][i])*INF + g->Link[i]->bw) ); //g->Link[i]->bw : residual bandwidth
		goal += cost[d]*req[d].flow;	
	}
	model.add(IloMaximize(env,goal));

	//流量守恒约束
	for(int d = 0; d < num; d++)
		for(int i = 0; i < g->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < g->adjL[i].size(); k++) // 出度边
				constraint += x[d][g->adjL[i][k]->id];
			for(unsigned int k = 0; k < g->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][g->adjRL[i][k]->id];
			// 出 - 入
			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}

		ORsolver.setOut(env.getNullStream());
		double obj = SMALL;
		ORsolver.solve();
		if(ORsolver.getStatus() == IloAlgorithm::Infeasible)
			env.out() << "throughput unfeasible" << endl;
		else{
			obj = ORsolver.getObjValue();

			// bw
			for(int i = 0; i < g->m; i++){  
				double load = 0;
				for(int d = 0; d < num; d++){
					load += ORsolver.getValue(x[d][i])*req[d].flow;
				}
				g->Link[i]->use = load;
			}
		}
		for(size_t i = 0; i < req.size(); i++)
			x[i].end();
		x.end();
		env.end();
		return obj;
}


double CGraph::EE(int id,int s,int t,double dm,bool needpath,double OPEN){
	vector<int> p, flag;
	vector<double> d;//energy
	for(int i = 0; i < n; i++){
		p.push_back(-1);
		flag.push_back(0);
		d.push_back(INF);
	}

	d[s] = 0;
	int cur = s;
	do{	
		flag[cur] = 1;
		for(unsigned int i = 0; i < adjL[cur].size(); i++){
			CEdge *e = adjL[cur][i];
			//e->latency = pow((e->use+dm),2)-pow(e->use,2);
			e->latency = pow((e->use+dm),2);
			if( 0 == e->use )
				e->latency += OPEN;
			if(e->capacity - e->use >= dm && d[e->head] > d[e->tail] + e->latency){
				d[e->head] = d[e->tail] + e->latency;
				p[e->head] = e->id;
			}
		}
		cur = -1;
		for(int i = 0; i < n; i++)
			if(!flag[i] && (cur == -1 || d[cur] > d[i] ))
				cur = i;
	}while(cur != -1);

	cur = t;
	do{
		if(p[cur] == -1)
			break;
		Link[p[cur]]->use += dm;
		cur = Link[p[cur]]->tail;
	}while(cur != s);

	if(needpath){
		reqPathID[id].clear();
		cur = t;
		do{
			if(p[cur] == -1)
				break;
			this->reqPathID[id].push_back(p[cur]);
			cur = Link[p[cur]]->tail;
		}while(cur != s);
		reverse(reqPathID[id].begin(),reqPathID[id].end());
	}

	return d[t];
}

#endif
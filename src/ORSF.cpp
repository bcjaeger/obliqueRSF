#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector surv_est(NumericVector times,
                       NumericVector probs,
                       NumericVector newtimes){
  
  int ntimes = times.length();
  NumericVector newprobs(newtimes.length());
  
  newprobs[newtimes<times[0]]=1;
  
  for(int i = 1; i<ntimes; i++){
    newprobs[(newtimes<times[i]) & (newtimes>=times[i-1])]=probs[i-1];
  }
  
  newprobs[newtimes>=times[ntimes-1]]=probs[ntimes-1];
  
  return(newprobs);
  
}


// [[Rcpp::export]]
double innerprod(NumericVector x,
                 NumericVector y){
  return(sum(x*y));
}

//[[Rcpp::export]]
double lrtestC(NumericVector time,
               NumericVector status,
               NumericVector grp){
  
  int n = time.length();
  double Y=n;
  double Y1=sum(grp);
  
  int lwr=0, upr=0;
  
  for(int i=0; i<n; i++){
    if(status[i]==0){
      upr++;
    } else {
      break;
    }
  }
  
  IntegerVector indx = seq(lwr,upr);
  double d = sum(as<NumericVector>(status[indx]));
  double d1= innerprod(as<NumericVector>(status[indx]),
                       as<NumericVector>(grp[indx])  );
  
  double e1=Y1*d/Y;
  double e0=(Y-Y1)*d/Y;
  double o1=d1;
  double o0=d-d1;
  
  double V = (Y-Y1)*Y1*d*(Y-d) / (pow(Y,2)*(Y-1));
  Y -= indx.length();
  Y1 -= sum(as<NumericVector>(grp[indx]));
  
  lwr=upr+1;
  int counter=1;
  
  for( ; ; ){
    
    upr=lwr;
    while((status[upr]==0) & (upr<n-1)){
      upr++;
    }
    
    if(upr==n-1){
      if(status[upr]==0){
        break;
      } else {
        
        indx=seq(lwr,upr);
        
        d=sum(as<NumericVector>(status[indx]));
        d1=innerprod(as<NumericVector>(status[indx]),
                     as<NumericVector>(grp[indx])  );
        
        e1+=(Y1*d/Y);
        e0+=((Y-Y1)*d/Y);
        o1+=d1;
        o0+=d-d1;
        
        V += (Y-Y1)*Y1*d*(Y-d) / (pow(Y,2)*(Y-1));
        Y -= indx.length();
        Y1 -= sum(as<NumericVector>(grp[indx]));
        
        break;
      }
    }
    
    indx=seq(lwr,upr);
    
    d=sum(as<NumericVector>(status[indx]));
    d1=innerprod(as<NumericVector>(status[indx]),
                 as<NumericVector>(grp[indx])  );
    
    e1+=(Y1*d/Y);
    e0+=((Y-Y1)*d/Y);
    o1+=d1;
    o0+=d-d1;
    
    V += (Y-Y1)*Y1*d*(Y-d) / (pow(Y,2)*(Y-1));
    Y -= indx.length();
    Y1 -= sum(as<NumericVector>(grp[indx]));
    counter++;
    lwr=upr+1;
    
    if(Y==1) break;
    
  }
  
  return pow(o1-e1,2) / V;
  
}

// [[Rcpp::export]]
NumericVector moving_average(NumericVector a){
  NumericVector result(a.length());
  result[0]=a[0];
  for(int i = 1; i<a.length(); i++){
    result[i] = result[i-1] + (a[i]-result[i-1]) / (i+1);
  }
  return(result);
}


// [[Rcpp::export]]
NumericVector colmeans(NumericMatrix input_mat){
  
  int ncol = input_mat.ncol();
  int nrow = input_mat.nrow();
  
  NumericVector output_vec(ncol);
  
  for(int i=0; i < ncol; i++){
    double total=0;
    for(int j=0; j < nrow; j++){
      total = total+input_mat(j,i);
    }
    output_vec[i]=total/nrow;
  }
  
  return(output_vec);
  
}


// [[Rcpp::export]]
String pick_node(double wt,
                 double cut_pnt,
                 CharacterVector options){
  
  if(wt <= cut_pnt){
    return options[0];
  } else {
    return options[1];
  }
  
}

// [[Rcpp::export]]
NumericMatrix predict_orsf(List forest,
                           NumericMatrix newx,
                           NumericVector times){
  
  List newx_attr = newx.attr("dimnames");
  CharacterVector ftrs = newx_attr[1];
  int ntree = forest.length();
  
  NumericMatrix output(newx.nrow(), times.length());
  
  for(int indx=0; indx < newx.nrow(); indx++){
    
    // create a vector of inputs for observation with given indx
    NumericVector vec = newx(indx,_); vec.names() = ftrs;
    NumericMatrix mat(ntree, times.length());
    
    for(int i = 0; i < ntree; i++) {
      
      List ostree = forest[i];
      
      // New observations start at the root of the tree
      List nodes = ostree["nodes"];
      List current_node = nodes["R"];
      bool is_leaf = current_node["leaf"];
      
      while(!is_leaf){
        
        // Identify the children of the current node
        CharacterVector children = current_node["children"];
        
        // Identify variables & coefficients for linear combos.
        CharacterVector bvrs = current_node["bvrs"];
        NumericVector vec_bvrs = vec[bvrs];
        NumericVector nod_bwts = current_node["bwts"];
        
        // compute the current observation's linear combination of inputs
        
        double wt = innerprod(vec_bvrs, nod_bwts);
        double ct_pnt = current_node["cut_pnt"];
        
        // Determine the next node for the current observation
        String new_node = pick_node(wt,ct_pnt,children);
        
        // Identify this node in the list of nodes
        current_node = nodes[new_node];
        
        // Assess whether this node is a leaf
        is_leaf = current_node["leaf"];
        
      }
      
      NumericVector node_times=current_node["times"];
      NumericVector node_probs=current_node["probs"];
      
      NumericVector preds = surv_est(node_times,
                                     node_probs,
                                     times);
      
      //Rcpp::Rcout << preds << std::endl;
      
      
      mat(i,_) = preds;
      
    }
    
    output(indx,_) = colmeans(mat);
    
  }
  
  return(output);
  
}


// [[Rcpp::export]]
NumericMatrix predict_ost(List ostree,
                          NumericMatrix newx, 
                          NumericVector times){
  
  List newx_attr = newx.attr("dimnames");
  CharacterVector ftrs = newx_attr[1];
  
  NumericMatrix output(newx.nrow(), times.length());
  
  for(int indx=0; indx < newx.nrow(); indx++){
    
    // create a vector of inputs for observation with given indx
    NumericVector vec = newx(indx,_);
    vec.names() = ftrs;
    
    // New observations start at the root of the tree
    List nodes = ostree["nodes"];
    List current_node = nodes["R"];
    bool is_leaf = current_node["leaf"];
    
    while(!is_leaf){
      
      // Identify the children of the current node
      CharacterVector children = current_node["children"];
      
      // Identify variables & coefficients for linear combos.
      CharacterVector bvrs = current_node["bvrs"];
      NumericVector vec_bvrs = vec[bvrs];
      NumericVector nod_bwts = current_node["bwts"];
      
      // compute the current observation's linear combination of inputs
      
      double wt = innerprod(vec_bvrs, nod_bwts);
      double ct_pnt = current_node["cut_pnt"];
      
      // Determine the next node for the current observation
      String new_node = pick_node(wt,ct_pnt,children);
      
      // Identify this node in the list of nodes
      current_node = nodes[new_node];
      
      // Assess whether this node is a leaf
      is_leaf = current_node["leaf"];
      
    }
    
    NumericVector node_times=current_node["times"];
    NumericVector node_probs=current_node["probs"];
    NumericVector preds = surv_est(node_times, node_probs, times);
    
    output(indx,_) = preds;
    
  }
  
  return(output);
  
}


// [[Rcpp::export]]
NumericVector seql(long double from, 
                   long double to, 
                   long unsigned int length_out){
  arma::vec x = arma::linspace(from, to, length_out); 
  return Rcpp::NumericVector(x.begin(), x.end());
}


// [[Rcpp::export]]
StringVector leaf_nodes(List nodes){
  
  StringVector leaves(nodes.length());
  IntegerVector keep(nodes.length());
  keep.fill(-1);
  
  for(int i=0; i<nodes.length(); i++){
    List node= nodes[i];
    if(node["leaf"]){
      String leaf_name = node["name"]; 
      leaves[i]=leaf_name;
      keep[i]=i;
    }
  }
  keep = keep[keep>=0];
  return leaves[keep];
}


// [[Rcpp::export]]
int unique_len_fast(NumericVector x){
  for(int i=1; i<x.length(); i++){
    if(x[i]!=x[0]){
      return(2);
    } 
  }
  return(1);
}

// [[Rcpp::export]]
IntegerVector filter_unique(IntegerVector v,
                            IntegerVector i){
  IntegerVector tmp = v[i];
  return unique(tmp);
}

// [[Rcpp::export]]
StringVector modify_string(StringVector x, 
                           String newchar,
                           IntegerVector newindx) {
  
  StringVector y=x;
  
  for(int i=0; i<newindx.length(); i++){
    int indx=newindx[i];
    y[indx]=newchar;
  }
  return y;
  
}

// [[Rcpp::export]]
StringVector namegen(String origin_name, StringVector rhs)
{
  StringVector lhs(rhs.length()); 
  lhs.fill(origin_name);
  R_xlen_t i = 0, sz = lhs.size();
  StringVector res(sz);
  
  for (std::ostringstream oss; i < sz; i++, oss.str("")) {
    oss << lhs[i] << rhs[i];
    res[i] = oss.str();
  }
  return res;
}

// [[Rcpp::export]]

IntegerVector find_nonconst_cols(NumericMatrix mat, 
                                 IntegerVector indx,
                                 int ncol){
  
  IntegerVector nonconst_cols(ncol);
  nonconst_cols.fill(-1);
  
  for(int i=0; i<ncol; i++){
    
    NumericVector mat_col = mat(_,i);
    NumericVector mat_vec = mat_col[indx];
    int nvals=unique_len_fast(mat_vec);
    if(nvals>1) nonconst_cols[i]=i;
  }
  
  return(nonconst_cols[nonconst_cols>=0]);
  
}

// [[Rcpp::export]]

IntegerVector soft_sample(IntegerVector x, int size){
  
  if(size > x.length()){
    return(x);
  } else {
    return(sample(x,size,false));
  }
  
}

// [[Rcpp::export]]
IntegerVector find_indx(StringVector node_ids,
                        String current_node){
  
  IntegerVector indx(node_ids.length());
  indx.fill(-1);
  
  for(int i=0; i<node_ids.length(); i++){
    if(node_ids[i]==current_node) indx[i]=i;
  }
  return(indx[indx>=0]);
}

// [[Rcpp::export]]

NumericVector comp_preds(IntegerVector& indx, 
                         int& node_nobs,
                         IntegerVector& node_cols,
                         NumericVector& bwts,
                         NumericMatrix& dmat){
  
  NumericVector preds(node_nobs);
  for(int k=0; k<node_nobs;k++){
    int r=indx[k];
    for(int j=0; j<node_cols.length(); j++){
      int c=node_cols[j];
      preds[k] += dmat(r,c)*bwts[j];
    }
  }
  
  return preds;
  
}

// [[Rcpp::export]]
List srv_R(NumericVector time_indx,
           IntegerVector status_indx,
           Function f){
  
  return as<List>(f(time_indx,status_indx));
  
}

// [[Rcpp::export]]
double lrt_R(NumericVector tmp_time,
             NumericVector tmp_grp,
             NumericVector tmp_status,
             Function f){
  
  return as<double>(f(tmp_time, tmp_grp, tmp_status));
  
}

// [[Rcpp::export]]
List boot_R(NumericMatrix dmat,
            NumericVector time,
            IntegerVector status,
            IntegerVector inb,
            Function f){
  
  return as<List>(f(dmat,time,status,inb));
  
}

// [[Rcpp::export]]
List eval_forest_R(NumericMatrix oob_predictions,
                   NumericVector time,
                   IntegerVector status,
                   NumericVector eval_times,
                   Function f){
  
  return as<List>(f(oob_predictions, time, status, eval_times));
  
}

// [[Rcpp::export]]
NumericMatrix net_R(Rcpp::NumericMatrix dmat,
                    Rcpp::NumericVector time,
                    Rcpp::IntegerVector status,
                    Rcpp::IntegerVector indx,
                    Rcpp::IntegerVector cols,
                    int dfmax,
                    NumericVector alpha,
                    Function f){
  
  return Rcpp::as<NumericMatrix>(f(dmat,time,status,indx,cols,dfmax,alpha));
  
}


// [[Rcpp::export]]
List tune_node(NumericMatrix& dmat, 
               NumericMatrix& bwts_mat,
               IntegerVector& indx,
               IntegerVector& node_cols,
               NumericVector& status_indx,
               NumericVector& time_indx,
               List& node, 
               int& node_nobs,
               int& min_obs_in_leaf_node,
               int& min_events_in_leaf_node,
               int& nsplit,
               double mincriterion){
  
  double best_lrstat=0;
  
  List output = List::create(
    Named("not_splittable")=true,
    _["bwts"]=NULL,
    _["cut_pnt"]=NULL,
    _["grp"]=NULL);
  
  for(int u=0; u<bwts_mat.cols(); u++){
    
    NumericVector tmp_bwts = bwts_mat(_,u);
    
    if(sum(abs(tmp_bwts))>0){
      
      NumericVector preds=comp_preds(indx, 
                                     node_nobs, 
                                     node_cols, 
                                     tmp_bwts, 
                                     dmat);
      
      NumericVector all_cut_pnts = unique(preds); 
      
      int ncuts = all_cut_pnts.length();
      int ncuts_to_sample;
      
      if(ncuts<=nsplit){
        ncuts_to_sample=ncuts;
      } else {
        ncuts_to_sample=nsplit;
      }
      
      NumericVector cut_pnts;
      
      if(ncuts>1){
        
        if(ncuts==2){
          cut_pnts=median(all_cut_pnts);
        } else {
          cut_pnts=sample(all_cut_pnts,ncuts_to_sample,false);
        }
        
        // Rcout<<std::endl<<"Cut points for "<<u+1<<
        //   " variable model:"<<cut_pnts<<std::endl;
        
        NumericVector::iterator cut_pnt;
        
        for(cut_pnt=cut_pnts.begin(); cut_pnt!=cut_pnts.end(); ++cut_pnt){
          
          double lrstat=0;
          
          //Rcout<<"Current cut: "<<*cut_pnt<<std::endl;
          NumericVector tmp_grp(preds.length()); tmp_grp.fill(1);
          LogicalVector set_to_zero=preds<=*cut_pnt;
          tmp_grp[set_to_zero]=0;
          
          int p1 = sum(tmp_grp);
          int p0 = tmp_grp.length() - sum(tmp_grp);
          
          NumericVector status_p0 = status_indx[tmp_grp==0];
          NumericVector status_p1 = status_indx[tmp_grp==1];
          
          double s0 = sum(status_p0);
          double s1 = sum(status_p1);
          
          if((p1>=min_obs_in_leaf_node) & 
             (p0>=min_obs_in_leaf_node) & 
             (s1>=min_events_in_leaf_node) & 
             (s0>=min_events_in_leaf_node)){
            
            lrstat=lrtestC(time_indx, status_indx, tmp_grp);
            
            //Rcout<<"LR stat for this cut: "<<lrstat<<std::endl;
            
            if(lrstat>mincriterion){
              output["not_splittable"]=false;
            }
            
            if(lrstat>best_lrstat){
              output["bwts"]=tmp_bwts;
              output["cut_pnt"]=*cut_pnt;
              output["grp"]=tmp_grp;
              best_lrstat=lrstat;
            }
            
          } 
          
        }
        
      }
      
    }
    
  }
  
  return output;
  
}

// [[Rcpp::export]]
List OST(NumericMatrix dmat,
         StringVector features,
         NumericVector alpha,
         NumericVector time,
         IntegerVector status,
         IntegerVector inbag_orsf_ids,
         int min_events_to_split_node,
         int min_obs_to_split_node,
         int min_obs_in_leaf_node,
         int min_events_in_leaf_node,
         int mtry,
         int dfmax,
         int nsplit,
         double mincriterion,
         Function surv_KM_Rfun,
         Function glmnet_Rfun){
  
  // initialize node_ids: a vector indicating
  // which node each observation belongs to
  
  int n = dmat.nrow();
  int p = dmat.ncol();
  
  StringVector node_ids(n); node_ids.fill("R");
  
  IntegerVector bootstrap_ids = unique(inbag_orsf_ids);
  
  List root = List::create(
    Named("name")="R",
    _["parent"]="root",
    _["depth"]=0,
    _["leaf"]=true,
    _["grow_me"]=true
  );
  
  List nodes = List::create(Named("R")=root);
  
  StringVector nodes_to_grow = leaf_nodes(nodes);
  StringVector nodes_in_queu = clone(nodes_to_grow);
  
  int stop_if_0 = 1;
  
  while(stop_if_0!=0){
    
    stop_if_0 = 0;
    StringVector::iterator ng;
    
    for(ng=nodes_to_grow.begin(); ng!=nodes_to_grow.end(); ++ng){
      
      String current_node_name = *ng;
      List current_node = nodes[current_node_name];
      
      // indx indicates which observations belong to the current node
      IntegerVector indx = find_indx(node_ids, current_node_name);
      int node_nobs = indx.length();
      
      StringVector tmp(1); tmp.fill(current_node_name);
      //Rcout<<"Growing "<<tmp<< ", n="<<node_nobs<<": ";
      
      // cols indicates columns with >1 unique value
      IntegerVector cols = find_nonconst_cols(dmat,indx,p);
      int node_ncols = cols.length();
      
      // event values for the current node
      IntegerVector status_indx = status[indx];
      NumericVector status_indx_nv = as<NumericVector>(status_indx);
      int node_nevents = sum(status_indx);
      
      // time values for the current node
      NumericVector time_indx = time[indx];
      
      // ids from the original dataset in the current node
      IntegerVector indx_ids = filter_unique(inbag_orsf_ids,indx);
      int node_nsubs = indx_ids.length();
      
      if(
        ((node_nsubs    >= min_obs_to_split_node) &
          (node_ncols   >= 1) &
          (node_nevents >= min_events_to_split_node) & 
          (Rcpp::as<bool>(current_node["grow_me"])) 
        )){
        
        stop_if_0++;
        
        IntegerVector node_cols = soft_sample(cols, mtry);
        
        StringVector bvrs = features[node_cols];
        
        // use glmnet here
        NumericMatrix bwts_mat=net_R(dmat,
                                     time,
                                     status,
                                     indx,
                                     node_cols,
                                     dfmax,
                                     alpha,
                                     glmnet_Rfun);
        
        
        List split_status = tune_node(dmat, 
                                      bwts_mat,
                                      indx,
                                      node_cols,
                                      status_indx_nv,
                                      time_indx,
                                      current_node, 
                                      node_nobs,
                                      min_obs_in_leaf_node,
                                      min_events_in_leaf_node,
                                      nsplit,
                                      mincriterion);
        
        if(Rcpp::as<bool>(split_status["not_splittable"])){
          
          current_node["grow_me"]=false;
          current_node["nsubs"]=node_nsubs;
          current_node["nevents"]=node_nevents;
          List node_srv = srv_R(time_indx,status_indx,surv_KM_Rfun);
          current_node["times"]=node_srv["times"];
          current_node["probs"]=node_srv["probs"];
          nodes[current_node_name]=current_node;
          nodes_in_queu = leaf_nodes(nodes);
          
        } else {
          
          
          current_node["bvrs"] = bvrs;
          current_node["bwts"] = split_status["bwts"];
          current_node["cut_pnt"] = split_status["cut_pnt"];
          
          StringVector child_names(2);
          child_names[0]="0";
          child_names[1]="1";
          
          StringVector new_names=namegen(current_node_name, child_names);
          NumericVector grp = as<NumericVector>(split_status["grp"]);
          
          for(int k=0; k<node_nobs; k++){
            int ik = indx[k];
            if(grp[k]==0){
              node_ids[ik]=new_names[0];
            } else {
              node_ids[ik]=new_names[1];
            }
            
          }
          
          IntegerVector node_table = table(node_ids);
          StringVector node_names = node_table.names();
          
          // Rcout<<"New nodes: ";
          // for(int v = 0; v<node_names.length();v++){
          //   if((node_names[v]==new_names[0]) | (node_names[v]==new_names[1])){
          //     Rcout << node_names[v] << " (n="  << node_table[v] << ") "; 
          //     
          //   }
          // }
          // Rcout<<std::endl;
          
          
          current_node["leaf"]=false;
          current_node["grow_me"]=false;
          current_node["children"] = new_names;
          int new_depth=current_node["depth"];
          new_depth++;
          
          String new_node_name0=new_names[0];
          String new_node_name1=new_names[1];
          
          List node0 = List::create(
            Named("name")=new_node_name0,
            _["parent"]=current_node_name,
            _["depth"]=new_depth,
            _["leaf"]=true,
            _["grow_me"]=true);
          
          List node1 = List::create(
            Named("name")=new_node_name1,
            _["parent"]=current_node_name,
            _["depth"]=new_depth,
            _["leaf"]=true,
            _["grow_me"]=true);
          
          nodes[new_node_name0]=node0;
          nodes[new_node_name1]=node1;
          nodes[current_node_name]=current_node;
          nodes_in_queu = leaf_nodes(nodes);
          
          
        }
        
        //nodes[current_node_name]=current_node;
        
      } else {
        
        // if the node does not qualify for further splitting
        
        // This can be used for print out diagnostics
        StringVector out_string(1);
        out_string.fill(current_node_name);
        //Rcout << "Not Grown: " << out_string << std::endl;
        
        current_node["grow_me"]=false;
        current_node["nsubs"]=node_nsubs;
        current_node["nevents"]=node_nevents;
        List node_srv = srv_R(time_indx,status_indx,surv_KM_Rfun);
        current_node["times"]=node_srv["times"];
        current_node["probs"]=node_srv["probs"];
        nodes[current_node_name]=current_node;
        nodes_in_queu = leaf_nodes(nodes);
        
        
        
      }
    }
    
    nodes_to_grow = nodes_in_queu;
    
  }
  
  
  List output = List::create(
    Named("nodes")=nodes,
    _["bootstrap_ids"]=bootstrap_ids
  );
  
  return(output);
  
}

//[[Rcpp::export]]
List ORSFcpp(NumericMatrix dmat,
             StringVector features,
             NumericVector alpha,
             NumericVector time,
             IntegerVector status,
             int min_events_to_split_node,
             int min_obs_to_split_node,
             int min_obs_in_leaf_node,
             int min_events_in_leaf_node,
             int mtry,
             int dfmax,
             int nsplit,
             int ntree,
             double mincriterion,
             bool verbose,
             Function surv_KM_Rfun,
             Function bootstrap_Rfun,
             Function glmnet_Rfun,
             Function forest_eval_Rfun){
  
  int n = dmat.nrow();
  IntegerVector orsf_ids=seq(0,n-1);
  
  IntegerVector inb;
  List forest(ntree), tree_dat(3);
  
  // Grow a forest
  
  for(int tree=0; tree < ntree; tree++){
    
    inb = sample(orsf_ids, n, true);
    
    tree_dat = boot_R(dmat,time,status,inb,bootstrap_Rfun);
    
    if(verbose){
      Rcout<<"Fitting tree no. "<< tree+1 << std::endl;
    }
    
    forest[tree]=OST(as<NumericMatrix>(tree_dat["mat"]),
                     features,
                     alpha,
                     as<NumericVector>(tree_dat["time"]),
                     as<IntegerVector>(tree_dat["status"]),
                     inb,
                     min_events_to_split_node,
                     min_obs_to_split_node,
                     min_obs_in_leaf_node,
                     min_events_in_leaf_node,
                     mtry, 
                     dfmax,
                     nsplit,
                     mincriterion,
                     surv_KM_Rfun,
                     glmnet_Rfun);
    
  }
  
  LogicalVector event_occurred = status==1;
  NumericVector event_times = time[event_occurred]; 
  long double   event_maxtime = max(event_times);
  NumericVector eval_times = seql(0,event_maxtime,10);
  NumericMatrix oob_pred(n, eval_times.length());
  
  // Compute out-of-bag predictions
  for(int row=0; row < dmat.nrow(); row++){
    
    // create a vector of inputs for observation with given row
    NumericVector vec = dmat(row,_);
    vec.names() = features;
    
    NumericVector oob_trees(ntree);
    oob_trees.fill(-1);
    int tree_counter = 0;
    
    for(int i=0; i<ntree; i++){
      List ostree = forest[i];
      NumericVector bstrap = ostree["bootstrap_ids"];
      bool in_bag = false;
      
      for(int j=0; j<bstrap.length(); j++){
        if(bstrap[j]==row){
          in_bag = true;
          j = bstrap.length();
        }
      }
      
      if(!in_bag){
        oob_trees[tree_counter]=i;
        tree_counter = tree_counter+1;
      }
    }
    
    oob_trees=oob_trees[oob_trees>=0];
    
    int ntree_oob = oob_trees.length();
    
    NumericMatrix mat(ntree_oob, eval_times.length());
    
    for(int i = 0; i < ntree_oob; i++) {
      
      int tree_row = oob_trees[i];
      List ostree = forest[tree_row];
      
      // New observations start at the root of the tree
      List nodes = ostree["nodes"];
      List current_node = nodes["R"];
      bool is_leaf = current_node["leaf"];
      
      while(!is_leaf){
        
        // Identify the children of the current node
        CharacterVector children = current_node["children"];
        
        // Identify variables & coefficients for linear combos.
        CharacterVector bvrs = current_node["bvrs"];
        NumericVector vec_bvrs = vec[bvrs];
        NumericVector nod_bwts = current_node["bwts"];
        
        // compute the current observation's linear combination of inputs
        
        double wt = innerprod(vec_bvrs, nod_bwts);
        double ct_pnt = current_node["cut_pnt"];
        
        // Determine the next node for the current observation
        String new_node = pick_node(wt,ct_pnt,children);
        
        // Identify this node in the list of nodes
        current_node = nodes[new_node];
        
        // Assess whether this node is a leaf
        is_leaf = current_node["leaf"];
        
      }
      
      NumericVector node_times=current_node["times"];
      NumericVector node_probs=current_node["probs"];
      
      NumericVector preds = surv_est(node_times,
                                     node_probs,
                                     eval_times);
      
      //Rcpp::Rcout << preds << std::endl;
      
      mat(i,_) = preds;
      
    }
    
    oob_pred(row,_) = colmeans(mat);
    
  }
  
  List oob_eval = forest_eval_Rfun(oob_pred,time,status,eval_times);
  
  List output = List::create(
    Named("forest")=forest,
    _["oob_error"]=oob_eval);
  
  return output;
  
  
}

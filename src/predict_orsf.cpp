#include <Rcpp.h>

using namespace Rcpp;
//' Survival estimate
//' 
//' this function is used to create survival step estimates.
//' 
//' @param times a numeric vector
//' @param probs a numeric vector
//' @param newtimes a numeric vector
//' @export
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

//' inner product
//' 
//' this function computes an inner product of two vectors.
//' 
//' @param y a numeric vector
//' @param x a numeric vector
//' @export
// [[Rcpp::export]]
double innerprod(NumericVector x,
                 NumericVector y){
  return(sum(x*y));
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

//' column means
//' 
//' this function computes column means of a matrix.
//' 
//' @param input_mat a numeric matrix
//' @export
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

//' pick a node
//' 
//' this function determines which node new observations go to
//' 
//' @param wt a numeric scalar
//' @param cut_pnt a numeric scalar
//' @param options a character vector of options
//' @export
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

//' which trees
//' 
//' this function identifies out-of-bag trees.
//' 
//' @param vec a numeric vector
//' @param forest a list of oblique survival trees
//' @export
// [[Rcpp::export]]
NumericVector which_trees(NumericVector vec,
                          List forest){
  
  int ntree = forest.length();
  
  bool internal = false;
  StringVector ftrs = vec.names();
  for(int i=0; i<ftrs.length(); i++){
    if (ftrs[i]=="orsf_id"){
      internal = true;
      i = ftrs.length();
    }
  }
  
  if(internal){
    
    double id = vec["orsf_id"];
    NumericVector which_Trees(ntree);
    int tree_counter = 0;
    
    for(int i=0; i<ntree; i++){
      List ostree = forest[i];
      NumericVector bstrap = ostree["bootstrap_ids"];
      bool in_bag = false;
      
      for(int j=0; j<bstrap.length(); j++){
        if(bstrap[j]==id){
          in_bag = true;
          j = bstrap.length();
        }
      }
      
      if(!in_bag){
        which_Trees[tree_counter]=i;
        tree_counter = tree_counter+1;
      }
      
    }
    
    NumericVector output(tree_counter);
    
    for(int i = 0; i<output.length(); i++){
      output[i]=which_Trees[i];
    }
    
    return(output);
    
  } else {
    
    NumericVector output(ntree);
    for(int i=0; i<ntree; i++){
      output[i]=i;
    }
    
    return(output);
    
  }

}

//' Predictions from an ORSF
//' @param forest A list of oblique survival trees.
//' @param newx a data matrix.
//' @param times A vector of times in the range of the response variable, e.g. times when the response is a survival object, at which to return the survival probabilities.
//' @return A matrix of survival probabilities containing 1 row for each observation and 1 column for each value in times.
//' @export
// [[Rcpp::export]]
NumericMatrix predict_orsf(List forest,
                           NumericMatrix newx,
                           NumericVector times){

  List newx_attr = newx.attr("dimnames");
  CharacterVector ftrs = newx_attr[1];
  
  NumericMatrix output(newx.nrow(), times.length());
  
  for(int indx=0; indx < newx.nrow(); indx++){
    
    // create a vector of inputs for observation with given indx
    NumericVector vec = newx(indx,_);
    vec.names() = ftrs;
    
    NumericVector oob_trees = which_trees(vec, forest);
    
    int ntree = oob_trees.length();
    
    NumericMatrix mat(ntree, times.length());
    
    for(int i = 0; i < ntree; i++) {
      
      int tree_indx = oob_trees[i];
      List ostree = forest[tree_indx];
      
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
      
      mat(i,_) = preds;
      
    }
    
    output(indx,_) = colmeans(mat);
    
  }
  
  return(output);
  
}

//' Predictions from an OST
//' @param ostree an oblique survival tree object.
//' @param newx a data matrix.
//' @param times A vector of times in the range of the response variable, e.g. times when the response is a survival object, at which to return the survival probabilities.
//' @return A matrix of survival probabilities containing 1 row for each observation and 1 column for each value in times.
//' @export
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



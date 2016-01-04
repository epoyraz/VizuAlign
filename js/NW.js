//String repeat function
String.prototype.str_out_times = function(max)   {
  if(max>1) {
    var ret = this;
    for(var i=0; i<max-1; i++){
        ret += this;
    }
  } else if(max==1){
    ret=''+this;
  } else if(max==0){
    ret='';
  }
  return ret;
}
//Get the score
var Get_score = function(SeqA_aligned,SeqB_aligned,Par_match,Par_mismatch,Par_gap,Par_gap_ext) {
	score=0;
	flag=0;
	for(i=0;i<SeqA_aligned.length;i++) {
		if(SeqA_aligned[i]===SeqB_aligned[i]) {
			score+=Par_match;
			flag=0;
		} else if(SeqA_aligned[i]==='_' || SeqB_aligned[i]==='_') {
			if(flag==0) {
				score-=Par_gap;
				flag=1;
			} else {
				score-=Par_gap_ext;
			}
		} else {
			score+=Par_mismatch;
			flag=0;
		}
	}
	return(score);
}
//Similarity function
var s = function(a,b,Par_match,Par_mismatch) {
	if(a===b) {
		return Par_match;
	} else {
		return Par_mismatch;
	}
}
//Initiate the linear matrix(A as row)
var linear_Get_Mat = function(SeqA,SeqB,d) {
	var mat=new Array();
	N=SeqA.length+1;
	M=SeqB.length+1;
  for(i=0;i<N;i++) {
		thisthing=[0];
    for(j=1;j<M;j++) {
      thisthing.push(0);
    }
    mat.push(thisthing);
	}
	for(j=0;j<M;j++) {
		mat[0][j]=0-j*d;
	}
	for(i=0;i<N;i++) {
		mat[i][0]=0-i*d;
	}
	return(mat);
}
//Linear trace back
var Linear_trace_back = function(SeqA,SeqB,Dir) {
	SeqA_aligned='';
	SeqB_aligned='';
	i=SeqA.length-1;
	j=SeqB.length-1;
	while (i>=0 && j>=0) {
		if(Dir[i+1][j+1]==0) {
			SeqA_aligned+=SeqA[i];
			SeqB_aligned+=SeqB[j];
			i-=1;
			j-=1;
		} else if(Dir[i+1][j+1]==1) {
			SeqA_aligned+='_';
			SeqB_aligned+=SeqB[j];
			j-=1;
		} else if(Dir[i+1][j+1]==-1) {
			SeqA_aligned+=SeqA[i];
			SeqB_aligned+='_';
			i-=1;
		}
	}
  var Linear_trace_back_res=new Object();
  Linear_trace_back_res.SeqA_aligned=SeqA_aligned;
  Linear_trace_back_res.SeqB_aligned=SeqB_aligned;
  Linear_trace_back_res.i=i;
  Linear_trace_back_res.j=j;
	return Linear_trace_back_res;
}
//Linear alignment
var Linear_alignment = function(SeqA,SeqB,Par_match,Par_mismatch,Par_gap) {
	//Get the recursion matrix
	Mat=linear_Get_Mat(SeqA,SeqB,Par_gap);
	Dir=linear_Get_Mat(SeqA,SeqB,0);	//0:from up-left, 1: from left, -1: from up
	for(i=0;i<SeqA.length;i++) {
		for(j=0;j<SeqB.length;j++) {
      sscore=s(SeqA[i],SeqB[j],Par_match,Par_mismatch);
      Mat[i+1][j+1]=Math.max(Mat[i][j]+sscore,Mat[i+1][j]-Par_gap,Mat[i][j+1]-Par_gap);
      if(Mat[i+1][j+1]==Mat[i][j]+sscore) {
        Dir[i+1][j+1]=0;
      } else if(Mat[i+1][j+1]==Mat[i][j+1]-Par_gap) {
        Dir[i+1][j+1]=-1;
      } else if(Mat[i+1][j+1]==Mat[i+1][j]-Par_gap) {
        Dir[i+1][j+1]=1;
      }
		}
	}
	//Trace back
  Linear_trace_back_res=Linear_trace_back(SeqA,SeqB,Dir);
	SeqA_aligned=Linear_trace_back_res.SeqA_aligned;
	SeqB_aligned=Linear_trace_back_res.SeqB_aligned;
  i=Linear_trace_back_res.i;
  j=Linear_trace_back_res.j;
  //Reverse the sequences
  SeqA_aligned=SeqA_aligned.split("").reverse().join("");
	SeqB_aligned=SeqB_aligned.split("").reverse().join("");
  if(i>=0) {
		SeqA_aligned=SeqA.slice(0,i+1)+SeqA_aligned;
		SeqB_aligned='_'.str_out_times(i+1)+SeqB_aligned;
	}
	if(j>=0) {
		SeqA_aligned='_'.str_out_times(j+1)+SeqA_aligned;
		SeqB_aligned=SeqB.slice(0,j+1)+SeqB_aligned;
	}
	//Get the score
	Score=Get_score(SeqA_aligned,SeqB_aligned,Par_match,Par_mismatch,Par_gap,Par_gap);
  var res = new Object();
  res.Mat = Mat;
  res.Dir=Dir;
  res.Score=Score;
  res.SeqA_aligned=SeqA_aligned;
  res.SeqB_aligned=SeqB_aligned;
  return res;
}
//Initiate the affine matrix(A as row)
var Get_Mat = function(SeqA,SeqB,d,e,Type) {
	mat=new Array();
	N=SeqA.length+1;
	M=SeqB.length+1;
	for(i=0;i<N;i++) {
		thisthing=[0];
    for(j=1;j<M;j++) {
      thisthing.push(0);
    }
    mat.push(thisthing);
	}
	if(Type==='M') {
		for(j=1;j<M;j++) {
			mat[0][j]=-999;
		}
		for(i=1;i<N;i++) {
			mat[i][0]=-999;
		}
	} else if(Type=='X') {
		for(j=0;j<M;j++) {
			mat[0][j]=-999;
		}
		for(i=1;i<N;i++) {
			if(i==1) {
				mat[i][0]=0-d;
			} else {
				mat[i][0]=mat[i-1][0]-e;
			}
		}
	} else if(Type==='Y') {
		for(j=1;j<M;j++) {
			if(j==1) {
				mat[0][j]=0-d;
			} else {
				mat[0][j]=mat[0][j-1]-e;
			}
		}
		for(i=0;i<N;i++) {
			mat[i][0]=-999;
		}
	} else if(Type==='D') {
		for(i=0;i<N;i++) {
			for(j=0;j<M;j++) {
				mat[i][j]=[0,0,0];
			}
		}
	}
	return(mat);
}
//Affine trace back
var Trace_back = function(SeqA,SeqB,SeqA_aligned,SeqB_aligned,Dir,thisthing,i,j) {
  while(i>=1 && j>=1) {
    newthing=Dir[i][j][thisthing];
		if(thisthing==0) {
      SeqA_aligned+=SeqA[i-1];
			SeqB_aligned+=SeqB[j-1];
		} else if(thisthing==2) {
      SeqA_aligned+='_';
			SeqB_aligned+=SeqB[j-1];
		} else if(thisthing==1) {
      SeqA_aligned+=SeqA[i-1];
			SeqB_aligned+='_';
		}
    if(thisthing==0) {
      i=i-1;
      j=j-1;
    } else if(thisthing==1) {
      i=i-1;
    } else if(thisthing==2) {
      j=j-1;
    }
    thisthing=newthing;
	}
  var Trace_back_res=new Object();
  Trace_back_res.SeqA_aligned=SeqA_aligned;
  Trace_back_res.SeqB_aligned=SeqB_aligned;
  Trace_back_res.i=i;
  Trace_back_res.j=j;
	return Trace_back_res;
}
//Affine alignment
var Alignment = function(SeqA,SeqB,Par_match,Par_mismatch,Par_gap,Par_gap_ext) {
	///Get the recursion matrix
	Mat=Get_Mat(SeqA,SeqB,Par_gap,Par_gap_ext,"M");
	Ix=Get_Mat(SeqA,SeqB,Par_gap,Par_gap_ext,"X");
	Iy=Get_Mat(SeqA,SeqB,Par_gap,Par_gap_ext,"Y");
	Dir=Get_Mat(SeqA,SeqB,1,1,"D");	//0:from up-left, 1: from left, -1: from up
  for(i=0;i<SeqA.length;i++) {
		for(j=0;j<SeqB.length;j++) {
			Dir[i+1][j+1]=[0,0,0];
			Ix[i+1][j+1]=Math.max(Mat[i][j+1]-Par_gap,Ix[i][j+1]-Par_gap_ext);
			if(Ix[i+1][j+1]==Mat[i][j+1]-Par_gap) {
				Dir[i+1][j+1][1]=0;
			} else {
				Dir[i+1][j+1][1]=1;
			}
			Iy[i+1][j+1]=Math.max(Mat[i+1][j]-Par_gap,Iy[i+1][j]-Par_gap_ext);
			if(Iy[i+1][j+1]==Mat[i+1][j]-Par_gap) {
				Dir[i+1][j+1][2]=0;
			} else {
				Dir[i+1][j+1][2]=2;
			}
			diag=Mat[i][j]+s(SeqA[i],SeqB[j],Par_match,Par_mismatch);
			left=Ix[i][j]+s(SeqA[i],SeqB[j],Par_match,Par_mismatch);
			up=Iy[i][j]+s(SeqA[i],SeqB[j],Par_match,Par_mismatch);
			Mat[i+1][j+1]=Math.max(diag,left,up);
			if(Mat[i+1][j+1]==diag) {
				Dir[i+1][j+1][0]=0;
			} else if(Mat[i+1][j+1]==left) {
				Dir[i+1][j+1][0]=1;
			} else if(Mat[i+1][j+1]==up) {
				Dir[i+1][j+1][0]=2;
			}
		}
	}
	//Start position of trace back
	Max=Math.max(Mat[SeqA.length][SeqB.length],Ix[SeqA.length][SeqB.length],Iy[SeqA.length][SeqB.length]);
	Fromlist=[];
	if(Max==Mat[SeqA.length][SeqB.length]) {
    Fromlist.push(0);
	}
	if(Max==Ix[SeqA.length][SeqB.length]) {
		Fromlist.push(1);
	}
	if(Max==Iy[SeqA.length][SeqB.length]) {
		Fromlist.push(2);
	}
  console.log(Fromlist)
  //Result list
  reslist=new Array();
	//Trace back
	for(m=0;m<Fromlist.length;m++) {
    var From;
    var SeqA_aligned;
    var SeqB_aligned;
    var i;
    var j;
    From=Fromlist[m];
    SeqA_aligned='';
    SeqB_aligned='';
    i=SeqA.length;
		j=SeqB.length;
    Trace_back_res=Trace_back(SeqA,SeqB,SeqA_aligned,SeqB_aligned,Dir,From,i,j);
    SeqA_aligned=Trace_back_res.SeqA_aligned;
    SeqB_aligned=Trace_back_res.SeqB_aligned;
    i=Trace_back_res.i-1;
    j=Trace_back_res.j-1;
		SeqA_aligned=SeqA_aligned.split("").reverse().join("");
		SeqB_aligned=SeqB_aligned.split("").reverse().join("");
    if(i>=0) {
			SeqA_aligned=SeqA.slice(0,i+1)+SeqA_aligned;
			SeqB_aligned='_'.str_out_times(i+1)+SeqB_aligned;
		}
		if(j>=0) {
			SeqA_aligned='_'.str_out_times(j+1)+SeqA_aligned;
			SeqB_aligned=SeqB.slice(0,j+1)+SeqB_aligned;
		}
    var Score
    Score=Get_score(SeqA_aligned,SeqB_aligned,Par_match,Par_mismatch,Par_gap,Par_gap_ext);
    var res = new Object();
    res.Mat = Mat;
    res.Ix=Ix;
    res.Iy=Iy;
    res.Dir=Dir;
    res.score=Score;
    res.SeqA_aligned=SeqA_aligned;
    res.SeqB_aligned=SeqB_aligned;
    reslist.push(res);
  }
  return reslist;
}

//affine: Alignment(SeqA,SeqB,Par_match,Par_mismatch,Par_gap,Par_gap_ext)
//linear: Linear_alignment(SeqA,SeqB,Par_match,Par_mismatch,Par_gap)

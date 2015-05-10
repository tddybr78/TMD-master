(* ::Package:: *)

(* ::Title:: *)
(*Mathematica package to manipulate and calculate with PDF files*)


(* ::Text:: *)
(*Version History :*)
(*6 Aug 2014-- Ben Clark, Creation of package to facilitate the use of LHA files with error calculation.*)
(*22 Sept 2014-- Eric Godat, updated usage and further combined LHAall into.*)
(*  *)


(* ::Section:: *)
(*Setup Package*)


BeginPackage["pdfCalc`"];


Print[" "];
Print[" - Required Package: pdfCalc --Loaded - "];
Print[" "];


(* ::Code:: *)
(***********************************************************************)

pdfFunction::usage="pdfFunction[setNumber,flavor,x,Q]: This function returns the value of the pdf for the .pds/.dat file \!\(\*
StyleBox[\"setNumber\",\nFontSlant->\"Italic\"]\), for 
the given flavor and value of Bjorken x and Q.

\!\(\*
StyleBox[\"Warning\",\nFontSlant->\"Italic\"]\): The results of this function are only reliable between the maximmum and minimum values of x and Q in the .pds/.dat file.";

(***********************************************************************)

pdfFunctionLHA::usage="pdfFunction[setNumber,flavor,x,Q]: This function returns the value of the pdf for the .pds/.dat file \!\(\*
StyleBox[\"setNumber\",\nFontSlant->\"Italic\"]\), for 
the given flavor and value of Bjorken x and Q.

\!\(\*
StyleBox[\"Warning\",\nFontSlant->\"Italic\"]\): The results of this function are only reliable between the maximmum and minimum values of x and Q in the .pds/.dat file.";

(***********************************************************************)


pdfLowFunction::usage="pdfLowFunction[setNumber,flavor,x,Q,[power]]: This function returns the value of the pdf as in pdfFunction, 
but with an exrapolation below the minimum x value that goes as \!\(\*FractionBox[\(1\), SuperscriptBox[\(x\), \(power\)]]\). The optional input power has a default value of power=1.0";
(***********************************************************************)

pdfSetList::usage="pdfSetList: This global variable prints the list of the loaded .pds/.dat grid files, their \!\(\*
StyleBox[\"setNumber\",\nFontSlant->\"Italic\"]\)., maximal number of quark flavors, 
and number of valence flavors in each PDF set.";

(***********************************************************************)

pdfResetCTEQ::usage="pdfResetCTEQ[]: This funtion deletes all .pds/.dat files from memory and resets all the internal variables in the package.

\!\(\*
StyleBox[\"Note\",\nFontSlant->\"Italic\"]\): This function does not accept any inputs.";

(***********************************************************************)

pdfFlavor::usage="pdfFlavor[flavor]: This function returns a string with the name of the PDF flavor. 

Example:
pdfFlavor[0], pdfFlavor[1], pdfFlavor[2] will return \"gluon\", \"up\", \"down\" for the gluon, up quark, and down quark PDFs.";

(***********************************************************************)

pdfXmin::usage="pdfXmin[setNumber]: This function returns the minimum x value in the PDF grid \!\(\*
StyleBox[\"setNumber\",\nFontSlant->\"Italic\"]\).";

(***********************************************************************)

pdfGetXlist::usage="pdfGetXlist[setNumber]: This function returns the grid values in x from the PDF set \!\(\*
StyleBox[\"setNumber\",\nFontSlant->\"Italic\"]\).";

(***********************************************************************)

pdfGetQlist::usage="pdfGetQlist[setNumber]: This function returns the grid values in Q from the PDF set \!\(\*
StyleBox[\"setNumber\",\nFontSlant->\"Italic\"]\).";

(***********************************************************************)

pdfAlphaS::usage="pdfAlphaS[setNumber,Q]:This function returns the value of \!\(\*SubscriptBox[\(\[Alpha]\), \(S\)]\) at hard scattering energy \!\(\*
StyleBox[\"Q\",\nFontSlant->\"Italic\"]\) when this information is available in the .pds or .info file. 

\!\(\*
StyleBox[\"Warning\",\nFontSlant->\"Italic\"]\): This function will print a text meesage if the \!\(\*SubscriptBox[\(\[Alpha]\), \(S\)]\) information is not available or if the file gives a single value for \!\(\*SubscriptBox[\(\[Alpha]\), \(S\)]\). 
The function then returns a Null value.";
(***********************************************************************)

pdfAlphaSLHA::usage="pdfAlphaSLHA[setNumber,Q]:This function returns the value of \!\(\*SubscriptBox[\(\[Alpha]\), \(S\)]\) at hard scattering energy \!\(\*
StyleBox[\"Q\",\nFontSlant->\"Italic\"]\) when this information is available in the .pds or .info file. 

\!\(\*
StyleBox[\"Warning\",\nFontSlant->\"Italic\"]\): This function will print a text meesage if the \!\(\*SubscriptBox[\(\[Alpha]\), \(S\)]\) information is not available or if the file gives a single value for \!\(\*SubscriptBox[\(\[Alpha]\), \(S\)]\). 
The function then returns a Null value.";
(***********************************************************************)




(* ::Section:: *)
(*define variables used in package*)


(* ::Input:: *)
(***********************************************************************)

nSetCount::usage="The number of PDF sets read into memory.";

(***********************************************************************)

pdfTableData::usage="Nested list containing the parsed data from the PDF sets.";

(***********************************************************************)

alphalist::usage="List containing the \!\(\*SubscriptBox[\(\[Alpha]\), \(S\)]\) data from the PDF sets.";

(***********************************************************************)

pdfxmin::usage="List of minimum values of Bjorken x for the PDF sets.";

(***********************************************************************)

(*TESTING THIS*)
numQpart::usage="number of partitions in a given set";

Begin["`Private`"]; 


(* ::Section:: *)
(*initialization*)


(* ::Input:: *)
(*(*this section will set up the global variables and should only be run once when this package is loaded. these variables are used by all packages.*)*)


(* ::Input:: *)
nSetCount=0;
pdfSetList={};
Clear[pdfTableData,alphalist,pdfxmin];
alphalist[_]=False;
pdfTableData[_]=False;(*this is for error control*)
pdfxmin[_]=False;
(* The above are global variables *)


(* ::Section:: *)
(*user control functions*)


(* ::Input:: *)
(*pdfResetCTEQ: no input, output is to reset memory and all global variables in the package*)
pdfResetCTEQ[]:=Module[{},
	nSetCount=0;
	pdfSetList={};
	Clear[pdfxmin,pdfTableData,alphalist];
alphalist[_]=False;
pdfxmin[_]=False;(*for error control*)pdfTableData[_]=False;(*for error control*)
	Print["All internal variables have been reset."]
  ];


(* ::Section:: *)
(*Interpolation*)


(* ::Subsection:: *)
(*This is a 4 - point (cubic) interpolation Function.*)



(*interpol4: a four-point Lagrange interpolation used throughout this package-->input is a point and a list (data) of the 4 points in x *)
(*           and Q surrounding the point,output is the value interpolated at the input point *)
interpol4[xIn_,data_]:=Module[
	{x,xVec,qVec,x0,x1,x2,x3,q0,q1,q2,q3,c0,c1,c2,c3,output},
	{xVec,qVec}=Transpose[data];
	x=xIn;
	{x0,x1,x2,x3}=xVec;
	{q0,q1,q2,q3}=qVec;
	c0=((x-x1)/(x0-x1)) ((x-x2)/(x0-x2)) ((x-x3)/(x0-x3));
	c1=((x-x0)/(x1-x0)) ((x-x2)/(x1-x2)) ((x-x3)/(x1-x3));
	c2=((x-x0)/(x2-x0)) ((x-x1)/(x2-x1)) ((x-x3)/(x2-x3));
	c3=((x-x0)/(x3-x0)) ((x-x1)/(x3-x1)) ((x-x2)/(x3-x2));
	output=q0 c0+q1 c1+q2 c2+q3 c3;
	If[output<0,Return[0],Return[output]]
  ];





(* ::Subsection:: *)
(*Interpolation Subfunctions*)



(*bigger: compares tow points xin and xintp and gives True if xin is greater than or equal to xintp*)
bigger[xin_,xintp_]:=(xin>=xintp);
SetAttributes[bigger,Listable];


(*findXindex:locates the point in a grid xlist that is immediately before the point x*)
(*input is point x and grid xlist, output is point in grid immdiately before x*)
findXindex[x_,xlist_]:=Module[{pos,xmin,xmax,r,xmis},
	xmin=xlist[[2]];
	r=Length[xlist];
	xmis=xlist[[r-1]];
	If[x<=xmin,Return[4],If[x>=xmis,Return[(r-1)],Return[Position[bigger[xlist,x],True]//First//First]];];
  ];


(*findQindex:locates the point in a grid qlist that is immediately before the point q*)
(*input is point x and grid qlist, output is point in grid immdiately before q*)
(*this is a unique function from findXindex as the structure of the two lists is different!!!*)
findQindex[q_,qlist_]:=Module[{pos,qmin,qmax,r,qmis},
	qmin=qlist[[2]];
	r=Length[qlist];
	qmis=qlist[[r-1]];
	If[q<=qmin,Return[3],If[q>=qmis,Return[(r-1)],Return[Position[bigger[qlist,q],True]//First//First]];];
  ];


(*getXdata && getQdata: construct the data lists for the interpolation for the appropriate points*)
(*these functions are unique*)
getXdata[xi_,qi_,flav_,xlist_,grid_]:=Transpose[{xlist[[xi-2;;xi+1]],grid[[flav,qi,xi-2;;xi+1]]}];
getQdata[xi_,qi_,flav_,qlist_,grid_]:=Transpose[{qlist[[qi-2;;qi+1]],grid[[flav,qi-2;;qi+1,xi]]}];


(*doQinterp: does the full interpolation over the 4 Q points starting from xi and qi for the flavor and grids *)
doQinterp[x_,xi_,qi_,flav_,xlist_,qlist_,grid_]:=
	Module[{pVec,qVec,output},
	pVec={
	interpol4[x,getXdata[xi,qi-2,flav,xlist,grid]],
	interpol4[x,getXdata[xi,qi-1,flav,xlist,grid]],
	interpol4[x,getXdata[xi,qi,flav,xlist,grid]],
	interpol4[x,getXdata[xi,qi+1,flav,xlist,grid]]};
	qVec=qlist[[qi-2;;qi+1]];
	output={qVec,pVec}//Transpose;
	Return[output]
  ];


(*fullinterp:does the interpolation over all the points in 4 directions outputs the value at that point*)
fullinterp[x_?NumericQ,q_?NumericQ,flav_,list_]:=
	Module[{xi,qi,qData,output,xlist,qlist,grid},
	xlist=list[[1]];
	qlist=list[[2]];
	grid=list[[3]];
	xi=findXindex[x,xlist];
	qi=findQindex[q,qlist];
	qData=doQinterp[x,xi,qi,flav,xlist,qlist,grid];
	output=interpol4[q,qData];
	Return[output];
  ];





(* ::Section:: *)
(*DefinePDF (CTEQ)*)


(* ::Subsection:: *)
(*The pdf function accepts values of x, q, parton #, and pdf data set in the following syntax: pdf[x,q,ipart,{dataset}].*)



(*plist: produced the value of the PDF at x, Q ,for flavor, and pdfDTabledata is the list*)
plist[x_?NumericQ,q_?NumericQ,flav_?IntegerQ,list_]:=fullinterp[x,q,flav,list];


(*pdfX: the internal version of the pdfFunction, input is set, flavor, x, and Q, output is the pdf*)
pdfX[iset_?IntegerQ,ipart_?IntegerQ,x_?NumericQ,q_?NumericQ]:=
	Module[{list,tpdf,nfmx,nfval,miscList,nfl},
	list=pdfTableData[iset];
	If[list==False,Return[Null]];  
	miscList=list[[4]];
	{nfl,nfmx,nfval}=miscList;
	If[Abs[ipart]>nfmx,Return[0]];
	If[ipart >= -nfmx && ipart <= nfval,
	  Return[plist[x,q,nfmx+ipart+1,list]],
	  Return[plist[x,q,nfmx-ipart+1,list]]
	];
  ];


pdfXLHA[iset_?IntegerQ,ipart_?IntegerQ,x_?NumericQ,q_?NumericQ]:=
	Module[{list,tpdf,nfmx,nfval,miscList,nfl,Qpartition},
	Qpartition=findQinpart[iset,q];
	list=tableDataMultiQ[iset,Qpartition];
	ipartmod=flavorshift[ipart];
	If[list==False,Return[Null]];  
	miscList=list[[4]];
	{nfl,nfmx,nfval}=miscList;
	Return[plist[x,q,ipartmod,list]];
  ];



(*pdfFunction: user version of the pdfX function deffinition*)
pdfFunction[iset_,ipart_,x_,q_] := pdfX[iset,ipart,x,q];

pdfFunctionLHA[iset_,ipart_,x_,q_] := pdfXLHA[iset,ipart,x,q];


(*pdfLowFunction: uset version of pdfX but with  an extrapolation below xmin that goes as x^-power*)
pdfLowFunction[iset_,ipart_,x_,q_,power_:1.0] := 
	Module[{xmin,output,list},
	If[pdfTableData[iset]==False,Return[Null]];  (*Check that table exists. 
	  Without this call function will not return null for a nonexistent list *)
	xmin=pdfxmin[iset];
	output=If[x>xmin,pdfX[iset,ipart,x,q],pdfX[iset,ipart,xmin,q] (xmin/x)^power];
	Return[output];
  ];





(* ::Section:: *)
(*Alpha_s implementation*)



(*getalphadata: constructs the data for the alphaS interoplation*)
getalphadata[qi_, alphalist_] := 
 Transpose[{alphalist[[1, qi - 2 ;; qi + 1]], 
   alphalist[[2, qi - 2 ;; qi + 1]]}
];


(*alphaS: input is set number and value of q, ouput is interpolated value of alphpaS at that point*)
alphaS[iset_,q_] := Module[{qini, alfa, data, output},
   alfa = Transpose[alphalist[iset]];
   qini = findQindex[q, alfa[[1]]];
   data = getalphadata[qini, alfa];
   output = interpol4[q, data];
   Return[output]
   ];

alphaSLHA[iset_,q_] := Module[{qini, alfa, data, output},
   alfa = alphalist[iset];
   qini = findQindex[q, alfa[[1]]];
   data = getalphadata[qini, alfa];
   output = interpol4[q, data];
   Return[output]
];

(*pdfAlphaSLHA: user available version of alphaS function. Will return various printed messages if it cannot interpolate*)
pdfAlphaSLHA[iset_,q_] := Module[{output},
    output = alphaSLHA[iset,q];
    Return[output];
];(*End pdfAlphaSLHA...*)


(*pdfAlphaS: user available version of alphaS function. Will return various prited messages if it cannot interpolate*)
pdfAlphaS[iset_,q_] := Module[{output},
   If[alphalist[iset] == False, 
    Return[Null]];(*Ensure that data exists*)  
   If[alphalist[iset] == {}, 
    Print["No \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\) information available for this .pds/.dat file."];
    Return[Null]];
   If[Length[alphalist[iset]] == 2, 
    Print["\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\) = ", alphalist[iset][[2]], " at q = ", alphalist[iset][[1]], 
     " GeV for all points in this .pds/.dat file."],
    output = alphaS[iset,q];
    Return[output]];
   ];






(*pdfFlavorTMP: a list (internal) of names of partons in CTEQ scheme*)
pdfFlavorTMP={
	"tbar",
	"bbar",
	"cbar",
	"sbar",
	"dbar",
	"ubar",
	"gluon",
	"up",
	"down",
	"strange",
	"charm",
	"bottom",
	"top"
	};


(*pdfFlavor:  User available function-->input is parton number in CTEQ scheme, output is string that is name of parton*)
pdfFlavor[i_]:=If[i!=21,pdfFlavorTMP[[i+7]],pdfFlavorTMP[[7]]
						];


(* ::Section:: *)
(*pdfGet Functions:*)



(*pdfGetXlist: input is set number, output is x grid from the associated PDF*)
pdfGetXlist[iset_]:=pdfTableData[iset][[1]] //Drop[#,1]&; (* The first x-point is for padding *)


(*pdfGetQlist: input is set number, output is Q grid from the associated PDF*)
pdfGetQlist[iset_]:=pdfTableData[iset][[2]];


(*pdfXmin: input is set number, ouput is the minimum value in the x grid in the assosiated PDF*)
pdfXmin[iset_] := Module[{value},
value=pdfxmin[iset];
If[value==False,Return[Null]];
Return[value];
];


(*pdfGetTableLHA: input is set number, output is grid values from the associated PDF*)
pdfGetTableLHA[iset_,Qpartition_:1]:=pdfTableDataMultiQ[iset,Qpartition][[3]];



(* ::Section:: *)
(*FUNCTIONS THAT ARE CLONES OF pdfParseLHAall FUNCTIONS.*)


(* ::Subsection:: *)
(*The pdf function accepts values of x, q, parton #, and pdf data set in the following syntax: pdf[x,q,ipart,{dataset}].*)



(*flavorshift: shifts the flavor numbers the user would input to refer to the correct flavor data in the grid*)
flist={1,2,3,5,4,11,7,6,8,9,10};(*The numbers refer to the column number in the .dat file and the order refers to the placement in numerical convention order -5\[Rule]5 with the gluon at 0*)
(*MODIFIED TO MATCH PDS IMPLEMENTATION 4<->5,6<->7*)
flavorshift[nf_]:=flist[[nf+6]];



(* ::Subsection:: *)
(*subfunctions*)


(* ::Code:: *)


(*Internal Function to adapt TableData for multiple q grids*)
tableDataMultiQ[nSetCount_,Qpartition_:1]:=Module[{tmp},
	tmp=Qpartition;
	If[tmp>numQpart[nSetCount],
		Print["Error: Q is not partitioned that many times in this set. 
				Consult the pdfNumQpartitionLHA function.\nShown is the value for the default value: Qpartition = 1"];
		 tmp=1
	];(*End If...*)
	xlist=pdfTableData[nSetCount][[1,tmp]];
	qlist=pdfTableData[nSetCount][[2,tmp]];
	list=pdfTableData[nSetCount][[3,tmp]];
	miscList=pdfTableData[nSetCount][[4,tmp]];
	Return[{xlist,qlist,list,miscList}];
];(*End pdfTabledataMultiQ...*)


(*Function that finds the number of q grids in a given data set*)
pdfNumQpartitionLHA[nSetCount_]:=numQpart[nSetCount];

(*Internal Function that takes a user inputed q value and locates which q partition that q falls into*)
findQinpart[iset_,q_]:=Module[{tmpqlist,numQ,i},
	numQ=numQpart[iset];
	If[numQ>1,Null,Return[1]];
	tmpqlist=pdfTableData[iset][[2]];
	For[i=1,
		i<=numQ,
		i++,
		If[q<=tmpqlist[[i,Length[tmpqlist[[i]]]]],
			If[q>= tmpqlist[[i,1]],
				Return[i];
				];(*End If...*)
			];(*End If...*)
		];(*End For...*)
	If[q<tmpqlist[[1,1]],
		Return[1],
		Return[numQ]
	];
];(*End pdfFindQ...*)




(* ::Section:: *)
(*End Package*)


(*Protect[pdfFunction,pdfLowFunction,pdfResetCTEQ,pdfFlavor,pdfXmin,pdfGetXlist,pdfGetQlist,
pdfAlphaS];*)(*protect all symbols*)

End[];  (* End Private Context *)

EndPackage[]; (* End Package Context *)

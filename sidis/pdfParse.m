(* ::Package:: *)

(* ::Title:: *)
(*Mathematica package to parse CTEQXX.pds files*)


(* ::Text:: *)
(*Original Version:*)
(*History : 15 Nov 2012-- v .0 .1, Ben Clark, Fred Olness*)
(*09 Jan 2013-- v .0 .2, Fred Olness, Olek Kusina Pass Iset, not whole list.Add miscList so we can deal with cases other than NF = 5.*)
(*15 Apr 2013-- v .0 .3, Pavel Nadolsky, implemented the new CT10 NNLO format, additional error control features, a flag for verbose printing in pdfParseCTEQ*)
(*20 Apr 2013-- v1 .0, Ben Clark, code reworked for consistancy and versions 0.1 and 0.2 merged, addtional error correction and alpha_s information included*)
(**)
(*New set of packages:*)
(*History: 30 October 2013- v1.1 Ben Clark, reworked to produce 3 packages. One to parse PDS files (This one), One to compute errors, and One to inport LHAgrid files for CTEQ*)
(*7 November 2013 v1.2 Ben Clark, rewrite whole package to do away with ifamily variable*)
(*6 Aug 2014 v1.3 Ben Clark, divide packages to produce a single package with all the common functions between the various packages.*)


(* ::Section:: *)
(*Setup Package*)


BeginPackage["pdfParse`","pdfCalc`"];(*Package requires pdfCalc.m for intrepolation routine. see http://ncteq.hepforge.org/code/pdf.html*)

Print["==============================================================="];
Print[" "];
Print[" - pdfParse - "];
Print["Version:  1.0"];
Print["Authors: D.B. Clark & F.I. Olness"];
Print[" "];
Print["Please cite: **************"];
Print["http://ncteq.hepforge.org/code/pdf.html"];
Print[" "];
Print["For a list of avialable commands, enter: ?pdf*"];
Print[" "];
Print["==============================================================="];


(***********************************************************************)

pdfFunction::usage="pdfFunction[setNumber,flavor,x,Q]: This function returns the value of the pdf for the .pds file setNumber, for 
the given flavor and value of Bjorken x and Q.

\!\(\*
StyleBox[\"Warning\",\nFontSlant->\"Italic\"]\): The results of this function are only reliable between the maximmum and minimum values of x and Q in the .pds file.";

(***********************************************************************)

pdfLowFunction::usage="pdfLowFunction[setNumber,flavor,x,Q,[power]]: This funtion returns the value of the pdf as in pdfFunction, 
but with an exrapolation below the minimum x value that goes as \!\(\*FractionBox[\(1\), SuperscriptBox[\(x\), \(power\)]]\). The optional input power has a default value of power=1.0";
(***********************************************************************)

pdfParseCTEQ::usage="pdfParseCTEQ[fileName,[verbose]]: This funciton reads an induvidual .pds file specified by fileName into memory. The 
optional parameter can be set to verbose=Flase to suppress the printout. The optional parameter has a default value verbose=True.";

(***********************************************************************)

pdfFamilyParseCTEQ::usage="pdfFamilyParseCTEQ[path,[fileType]]: This function reads all the files of type $fyleType$ in the 
directory path and stores them in memory. The function returns a list of set numbers that can be used to define a list. 
The optional input fileType has a default value of \"*.pds\".

Example:
  pdfFamilyParseCTEQ[\"MyGrids\",\"ct10*pds\"] reads all .pds files in the subdirectory \"MyGrids\" beginning with \"ct10\" into memory.";

(***********************************************************************)


Begin["`Private`"]; 


(* ::Section:: *)
(*Main Parsing Functions*)


(*typer: input is a PDS file, output is a value to control the behaviour of the parseing routine (ipdsformat).*)
typer[stream_] := Module[{chk0, chk1, ipdsformat},
   SetStreamPosition[stream, 0];
   Read[stream, String];
   chk0 = StringTake[Read[stream, Record, RecordSeparators -> ","], -3];
   Read[stream, Record];
   Read[stream, String];
   chk1 = StringTake[Read[stream, Record, RecordSeparators -> ","], -3];
   ipdsformat = Null; (* Use this if nothing matches *)
   If[chk0 == "rdr" && chk1 != "fMx",ipdsformat = 6];(*CTEQ6 .6.pds format;alpha_s is not specified*)
   If[chk0 == "ipk", ipdsformat = 12];(*CT12.pds format*)
   If[chk1 == "Ipk", ipdsformat = 10 ];(*Pre-CT10 NNLO format*)
   If[chk0 == "rdr" && chk1 == "fMx", ipdsformat = 99];(*LHA->PDS format for internal conversion routine*)
   Return[ipdsformat];
   ];



(* ::Code:: *)
(********************************** MAIN PARSING ROUTINE**********************************)

(*pdfParseCTEQ: input is a .pds file from cteq 6.6, 9, 10, 10nlo, 10nnlo, and 12 .pds files*)
(*output is the data in memory in an array called pdfTableData-->files added sequentialy and identified by nSetCount*)
pdfParseCTEQ[filename_?StringQ,verbose_:True]:=
	Module[{stream,order,ipk,nfl,qalpha,alfaQ,lambda,m1,m2,m3,m4,m5,m6,
	ipd0,ihdn,iknl,nfmx,nfval,nx,nt,ng,dum,idum,qini,qmax,xmin,xcr,xlist,nblk,
	length,biglist,nx1,nt1,list1,xlist1,qlist,nflav,qlist1,pdstype,tmpstring,
	qbase,aimass,fswitch,readlist,alphalist0,miscList,nn,end,bigtmp0,bigtmp,biglist0,str,endpoint},

	stream=OpenRead[filename];
	pdstype=typer[stream];
    SetStreamPosition[stream,0];
	tmpstring=Read[stream,String];(*name of PDF set*)
	Read[stream,String];

	If[pdstype==6,{order,nfl,lambda,m1,m2,m3,m4,m5,m6}=Read[stream,Table[Number,{i,1,9}]];
	  Read[stream,String];
  	{ipd0,ihdn,iknl,nfmx,nfval,dum,dum}=Read[stream,Table[Number,{i,1,7}]]
	];(*cteq 6*)
	If[pdstype==10,{order,nfl,qbase,m1,m2,m3,m4,m5,m6}=Read[stream,Table[Number,{i,1,9}]];
	  Read[stream,String];
	  {ipk,alfaQ,qalpha,nfmx,nfval,dum,dum}=Read[stream,Table[Number,{i,1,7}]]
	];(*cteq 10*)
	If[pdstype==12,{ipk,order,qalpha,alfaQ,m1,m2,m3,m4,m5,m6}=Read[stream,Table[Number,{i,1,10}]];
	  Read[stream,String];
	  {aimass,fswitch,ipd0,ihdn,iknl,nfmx,nfval}=Read[stream,Table[Number,{i,1,7}]]
	];(*cteq 12*)
	If[pdstype==99,{order,alfaQ,m1,m2,m3,m4,m5,m6}=Read[stream,Table[Number,{i,1,8}]];
	  Read[stream,String];
	  {nfmx,nfval}=Read[stream,{Number,Number}];
	  ng=0(*this value is not in the header*)
	];(*LHA->PDS format*)

	Read[stream,String];
	If[pdstype!=99,
	  {nx,nt,idum,ng,idum}=Read[stream,Table[Number,{i,1,5}]],
	  {nx,nt}=Read[stream,{Number,Number}]
	];

	If[ng>0, Read[stream,Table[Record,{i,1,ng+1}]]];
	Read[stream,Record];
    {qini,qmax}=Read[stream,{Number,Number}];

	If[pdstype==12||pdstype==99,(*Post-CT10 or LHA->PDS format *)
	  readlist=Read[stream,Table[LF[Number,Number,Number],{i,0,nt}]];
	  qlist=readlist/.LF[a__]:>{{a}[[1]],{a}[[2]]};
	  alphalist0=readlist/.LF[a__]:>{{a}[[1]],{a}[[3]]},
	  (*Pre-CT10 format *)
	  qlist=Read[stream,Table[{Number,Number},{i,0,nt}]]
	];

	Read[stream,String];
	If[pdstype!=99,{xmin,xcr}=Read[stream,{Number,Number}];
	  xlist=Read[stream,Table[Number,{i,1,nx}]],
	  xmin=Read[stream,Number];
	  xlist=Read[stream,Table[Number,{i,1,nx+1}]]
	];
	Read[stream,String];
	nblk=(nx+1)*(nt+1);
	length=nblk*(nfmx+1+nfval);
	biglist=Read[stream,Table[Number,{i,1,length}]];(*read in the PDF list*)
	endpoint=Read[stream];

	If[endpoint!=EndOfFile,
	  Print["Warning--Parton Distribution does not match header information for ",tmpstring,". Results are not reliable. 
			Please report file to the Developers"]
	];(*Check to see that you are at the end of the file*)

	nx1=nx+1;(*number of x elements*)
	nt1=nt+1;(*number of q elements*)
	nflav=(nfmx+1+nfval);(*number of flavors*)
	list1=Partition[Partition[biglist,{nx1}],{nt1}];(*partition the PDF file by x then q then flavor is left over*)
	If[pdstype!=99,xlist1=Join[{0.0},xlist],xlist1=xlist]; (*this adds a lower bound to the xlist*)
	qlist1=qlist//Transpose//First;
	Close[stream];

	If[pdstype==6,miscList={nfl,nfmx,nfval}]; (* CAREFUL: THE CONTENTS OF THIS MAY CHANGE!!! *)
	If[pdstype==10,miscList={nfl,nfmx,nfval}];
	If[pdstype==12||pdstype==99,miscList={nfl,nfmx,nfval}];


(**THIS ROUTINE IS WHERE WE WRITE THE LIST TO THE DATA TABLE**)
	If[Length[$MessageList] == 0,
    If[(verbose!=False),Print[tmpstring]];(*end if...*)(*verbose behavior of parsing routine*)

    nSetCount = nSetCount + 1;
    pdfxmin[nSetCount]=xmin;  (*internal variable--set so that funciton does not have to reference list each time*)
    pdfTableData[nSetCount]={xlist1,qlist1,list1,miscList};
    
    If[pdstype==6,alphalist[nSetCount]={}]; 
    If[pdstype==10,alphalist[nSetCount]={qalpha,alfaQ}];
    If[pdstype==12||pdstype==99,alphalist[nSetCount]=alphalist0];(*set up alphalist*)

    AppendTo[pdfSetList,{nSetCount,filename,nfmx,nfval}];
    Return[nSetCount], 
    Print[filename," was not initialized: ",Length[$MessageList]," error messages"];(*THIS IS THE ELSE STATEMENT*)
    
    Unprotect[$MessageList];
    Do[$MessageList = {}, {1}];
    Protect[$MessageList](*this will reset the $MessageList for the next file*)
	];(*If[Length[$MessageList]==0*)
  ];(*END OF PARSE FUNCTION*)

(********************************** END MAIN PARSING ROUTINE**********************************)






(*pdfFamilyParseCTEQ: input is a directory containing PDS files (optional input is filetype-->default is "*.pds")*)
(*                     all PDS files are stored in memory and output is a list containing the setNumbers of the PDS files that were stored *)
pdfFamilyParseCTEQ[path_?StringQ,type_:"*.pds"]:=
	Module[{currDir,fileList,output,fLength,startSet,endSet},
	currDir=Directory[];
	SetDirectory[path];
	fileList=FileNames[path<>"/"<>type];
	fLength=Length[fileList];
	startSet=nSetCount+1;
	Do[pdfParseCTEQ[fileList[[i]],False],{i,1,fLength}];
	endSet=nSetCount;
	Print["Included ",endSet-startSet+1," files in the PDF family."];
	SetDirectory[currDir];
	output=Range[startSet,endSet];
	Return[output];
  ];





(* ::Section:: *)
(*End Package*)


Protect[pdfFunction,pdfLowFunction,pdfParseCTEQ,pdfFamilyParseCTEQ];

End[];  (* End Private Context *)

EndPackage[]; (* End Package Context *)

/*
 * Result.cpp
 *
 *  Created on: Sep 26, 2013
 *      Author: jinzhang
 */

#include "Result.h"



Result::Result() {
	// TODO Auto-generated constructor stub

}

Result::~Result() {
	// TODO Auto-generated destructor stub
}

int Result::addResult(result_t result) {
	results.push_back(result);
	return 0;
}
//For troubleshooting
int Result::getResultSize() {
	cout<<"Results size:"<<results.size()<<endl;
	for(int k=0;k<results.size();k++)
		cout<<k<<" "<<results[k].nm5p<<" "<<results[k].nm3p<<endl;
	return 0;
}

int Result::searchResult(int geneId1, int geneId2, result_t& result) {
	for(int i=0;i<results.size();i++)
	{
		if((results[i].geneId1==geneId1 && results[i].geneId2==geneId2) || (results[i].geneId1==geneId2 && results[i].geneId2==geneId1))
		{
			result=results[i];
			return 1;
		}
	}
	return 0;
}


string getType(int tp)
{
    switch(tp)
    {
        case 0:
            return "Inter_Chromosomal_Fusion";
            break;
        case 1:
            return "Intra_Chromosomal_Fusion";
            break;
        case 2:
            return "Read_Through";
            break;
		case 3:
			return "Splice_Variant_Junction";
			break;
		case 4:
			return "Backsplice_Junction";
			break;
		case 5:
			return "Unknown_Junction";
			break;
/*       
 case 3:
            return "Overlap_Converging";
            break;
        case 4:
            return "Overlap_Diverging";
            break;
        case 5:
        	return "Adjacent_Converging";
        	break;
        case 6:
        	return "Adjacent_Diverging";
        	break;
*/
        default:
            return "Error";
            break;
    }

    return 0;
}

int printOneEncompassDna(encompass_dna_t &et, Reference & ref, ofstream & outFile)
{

	outFile<<et.name;
	outFile<<"\t";
	if(et.strand1==0)
		outFile<<"+";
	else
		outFile<<"-";
	outFile<<"\t";
	outFile<<ref.getCharName(et.tid1);
	outFile<<"\t";
	outFile<<et.pos1;
	outFile<<"\t";
	outFile<<et.len1;
	outFile<<"\t";
	if(et.strand2==0)
		outFile<<"+";
	else
		outFile<<"-";
	outFile<<"\t";
	outFile<<ref.getCharName(et.tid2);
	outFile<<"\t";
	outFile<<et.pos2;
	outFile<<"\t";
	outFile<<et.len2;
	outFile<<"\t";
	for(int i=0;i<et.seq1.size();i++)
		outFile<<et.seq1[i];
	outFile<<"\t";
	for(int i=0;i<et.seq2.size();i++)
		outFile<<et.seq2[i];
	outFile<<endl;
	return 0;
}


int printOneSplitDna(split_dna_t &st,Reference & ref,ofstream & outFile) {


	outFile<<st.name<<"\t";
	if(st.strand1==0)
		outFile<<"+\t";
	else
		outFile<<"-\t";
	outFile<<ref.getCharName(st.tid1)<<"\t";
	outFile<<st.pos1<<"\t";
	outFile<<st.len1<<"\t";

	if(st.strand2==0)
		outFile<<"+\t";
	else
		outFile<<"-\t";
	outFile<<ref.getCharName(st.tid2)<<"\t";
	outFile<<st.pos2<<"\t";
	outFile<<st.len2<<"\t";

	for(int j=0;j<st.seq.size();j++)
	{
		outFile<<st.seq[j];
	}
	outFile<<endl;

	return 0;

}

int outIndex;

int Result::printOneResult(int index, ofstream & outFile, Reference & ref, int isRunningNormal) {
	result_t rt=results[index];


	outFile<<"Fusion Candidate "<<++outIndex<<" : ";
	if(rt.isReci==1)
	{
		outFile<<rt.nm5p<<"<>"<<rt.nm3p;
	}
	else
	{
		outFile<<rt.nm5p<<">>"<<rt.nm3p;
	}

	outFile<<" NUM_EN_RNA "<<rt.numOfEnRna<<" NUM_SP_RNA "<<rt.numOfSpRna;
	if(indi > 1)
	{
		if(isRunningNormal==0)
			outFile<<" NUM_EN_DNA_Tumor "<<rt.numOfEnDnaT<<" NUM_SP_DNA_Tumor "<<rt.numOfSpDnaT;
		else
			outFile<<" NUM_EN_DNA_Normal "<<rt.numOfEnDnaT<<" NUM_SP_DNA_Normal "<<rt.numOfSpDnaT;

		if(indi>2)
		{
			if(isRunningNormal==0)
				outFile<<" NUM_EN_DNA_Normal "<<rt.numOfEnDnaN<<" NUM_EN_RNA_Normal "<<rt.numOfSpDnaN<<endl;
			else
				outFile<<" NUM_1 "<<rt.numOfEnDnaN<<" NUM_2 "<<rt.numOfSpDnaN<<endl;

		}
		else
		{
			outFile<<endl;
		}
	}
	else
	{
		outFile<<endl;
	}



	outFile<<"Encompassing RNA: "<<rt.enrnas.size()<<endl;
	for(int i=0;i<rt.enrnas.size();i++)
	{
		encompass_rna_t et=rt.enrnas[i];
		outFile<<et.name<<"\t";
		if(et.strand1==0)
			outFile<<"+\t";
		else
			outFile<<"-\t";
		outFile<<ref.getCharName(et.tid1)<<"\t";
		outFile<<et.pos1<<"\t";
		if(et.strand2==0)
			outFile<<"+\t";
		else
			outFile<<"-\t";
		outFile<<ref.getCharName(et.tid2)<<"\t";
		outFile<<et.pos2<<"\t";

		for(int j=0;j<et.seq1.size();j++)
		{
			outFile<<et.seq1[j];
		}
		outFile<<"\t";

		for(int j=0;j<et.seq2.size();j++)
		{
			outFile<<et.seq2[j];
		}
		outFile<<"\t";

		outFile<<et.numCopy;
		outFile<<endl;
	}

	outFile<<"Spanning RNA: "<<rt.numOfSpRna<<endl;

	int start=0;
	for(int i=0;i<rt.types.size();i++)
	{
		outFile<<"Splicing "<<i+1<<" : ";
		if(rt.primeOKs[i]==1)
		{
			outFile<<rt.nm5p<<">>"<<rt.nm3p;
		}
		else
		{
			outFile<<rt.nm3p<<">>"<<rt.nm5p;
		}
		outFile<<" ";
		outFile<<getType(rt.types[i])<<" ";
		if(rt.canos[i]==1)
			outFile<<"Canonical"<<" ";
		outFile<<rt.numOfsps[i]<<endl;

		for(int j=start;j<start+rt.numOfsps[i];j++)
		{
			outFile<<rt.sprnas[j].name<<"\t";
			if(rt.sprnas[j].strand1==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(rt.sprnas[j].tid1)<<"\t"<<rt.sprnas[j].pos1<<"\t"<<rt.sprnas[j].len1<<"\t";
	      //  {
	      //      outFile<<"*"<<" "<<"*"<<" "<<"*"<<" "<<"*"<<" ";
	      //  }
			if(rt.sprnas[j].strand2==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
	        outFile<<ref.getCharName(rt.sprnas[j].tid2)<<"\t"<<rt.sprnas[j].pos2<<"\t"<<rt.sprnas[j].len2<<"\t";
	        for(int x=0;x<rt.sprnas[j].seq.size();x++)
	        {
	        	outFile<<rt.sprnas[j].seq[x];
	        }
	        outFile<<"\t"<<rt.sprnas[j].hits<<endl;
		}

		start+=rt.numOfsps[i];
	}

	if(indi>1)
	{
		if(isRunningNormal==0)
			outFile<<"Encompassing DNA Tumor: "<<rt.endna1.size()<<endl;
		else
			outFile<<"Encompassing DNA Normal: "<<rt.endna1.size()<<endl;
		for(int i=0;i<rt.endna1.size();i++)
		{
			printOneEncompassDna(rt.endna1[i], ref, outFile);
		}
		if(isRunningNormal==0)
			outFile<<"Spanning DNA Tumor: "<<rt.spdna1.size()<<endl;
		else
			outFile<<"Spanning DNA Normal: "<<rt.spdna1.size()<<endl;

		for(int i=0;i<rt.spdna1.size();i++)
		{
			printOneSplitDna(rt.spdna1[i], ref, outFile);
		}

		if(indi>2)
		{
			if(isRunningNormal==0)
				outFile<<"Encompassing DNA Normal: "<<rt.endna2.size()<<endl;
			else
				outFile<<"Something Encompassing: "<<rt.endna2.size()<<endl;
			for(int i=0;i<rt.endna2.size();i++)
			{
				printOneEncompassDna(rt.endna2[i], ref, outFile);
			}
			if(isRunningNormal==0)
				outFile<<"Spanning DNA Normal: "<<rt.spdna2.size()<<endl;
			else
				outFile<<"Something Spanning: "<<rt.spdna2.size()<<endl;

			for(int i=0;i<rt.spdna2.size();i++)
			{
				printOneSplitDna(rt.spdna2[i], ref, outFile);
			}
		}

	}


	return 0;
}

//Goes to reads.txt by default
int Result::printAllResult(char* file, Reference & ref, int isRunningNormal, char* dir) {
	std::string filename = createFilename(dir, file);
	ofstream outFile(filename);
	outIndex=0;
	for(int i=0;i<results.size();i++)
	{
		if(results[i].realPrint==0)
			continue;
		printOneResult(i,outFile, ref, isRunningNormal);
	}
	return 0;
}




int Result::getTiers(double pn) {

	for(int i=0;i<results.size();i++)
	{
		result_t * prt=&(results[i]);

		if(indi==3)
		{
			int nn=(prt->numOfEnDnaN+prt->numOfSpDnaN);
			int nt=(prt->numOfEnDnaT+prt->numOfSpDnaT);
			if(nn>0)
			{
				int isNormalReal=0;
				if(nn>=nt)
				{
					isNormalReal=1;
				}
				else
				{
					if((double)nn/(double)nt>pn)
						isNormalReal=1;
				}
				if(isNormalReal==1)
				{
					prt->tier=7;
					continue;
				}

			}

		}

		//not 7;

		if(indi==1)
		{
			if(prt->isCanonical==1)
				prt->tier=3;
			else
				prt->tier=6;
			continue;
		}

		if(prt->isCanonical==1)//1 2 3
		{
			if(prt->numOfEnDnaT>0 && prt->numOfSpDnaT>0)
			{
				prt->tier=1;
				continue;
			}
			else if(prt->numOfEnDnaT>0)
			{
				prt->tier=2;
				continue;
			}
			else
			{
				prt->tier=3;
				continue;
			}

		}
		else// 4 5 6
		{
			if(prt->numOfEnDnaT>0 && prt->numOfSpDnaT>0)
			{
				prt->tier=4;
				continue;
			}
			else if(prt->numOfEnDnaT>0)
			{
				prt->tier=5;
				continue;
			}
			else
			{
				prt->tier=6;
				continue;
			}
		}
	}


	return 0;
}

result_t * Result::getOneResult(int index) {
	result_t * prt;
	if(index<0 || index>=results.size())
	{
		cerr<<"trying to find a result record that does not exist."<<endl;
		exit(0);
	}
	else
	{
		prt=&(results[index]);
	}
	return prt;
}

bool my_sort_result_func(result_t i, result_t j)
{
	if(i.tier<j.tier)
	{
		return true;
	}
	else if(i.tier==j.tier)
	{
		int t1=6;
		int t2=6;
		for(int x=0;x<i.types.size();x++)
		{
			if(i.types[x]<t1)
			{
				t1=i.types[x];
			}
		}
		for(int x=0;x<j.types.size();x++)
		{
			if(j.types[x]<t2)
			{
				t2=j.types[x];
			}
		}

		if(t1<t2)
		{
			return true;
		}
		else if(t1==t2)
		{
			if((i.numOfSpRna+i.numOfEnRna)>(j.numOfSpRna+j.numOfEnRna))
				return true;
			else if((i.numOfSpRna+i.numOfEnRna)==(j.numOfSpRna+j.numOfEnRna))
			{
				if(i.nm5p.compare(j.nm5p)<0)
				{
					return true;
				}
				else if(i.nm5p.compare(j.nm5p)==0)
				{
					if(i.nm3p.compare(j.nm3p)<0)
					{
						return true;
					}
					else
					{
						return false;
					}
				}
				{
					return false;
				}
			}
			else
				return false;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}

}


bool my_sort_int(int i, int j)
{
	return i<j;
}
bool my_sort_string(string i, string j)
{
	return i.compare(j);
}

int Result::removeMultiple(Gene & g, int largeNum)
{
	//default of largeNum == 4
	//If a gene shows up more than largeNum times, remove it
//cout<<"in removeMultiple"<<endl;

	vector<int> badGeneIds;
	vector<int> badGeneIds2;

	vector<int> badValues;
	vector<int> gids;

	for(int i=0;i<results.size();i++)
	{
		result_t * prt=this->getOneResult(i);
		gids.push_back(prt->geneId1);
		gids.push_back(prt->geneId2);
	}
	sort(gids.begin(),gids.end(),my_sort_int);
//cout<<"compute gids and value"<<endl;
//gids is a vector of all geneIds for all results
	int diffp=0;
	for(int i=1;i<gids.size();i++)
	{
		//If curr id differs from prev id, then
		//we are on a new id. diffp=1 means the 
		//curr id will not get added to badGeneIds
		if(gids[i]!=gids[i-1])
			diffp=i; 
		int end=0;
		if(i==gids.size()-1) //if last id, possible bad id
		{
			end=1;
		}
		else if(gids[i]!=gids[i+1]) //if new id, possible bad id
		{
			end=1;
		}
		//if last id or new id and
		if(end==1 && i-diffp+1>=largeNum)//para
		{
			badGeneIds.push_back(gids[i]);
			badValues.push_back(i-diffp+1);
		}
	}
//cout<<"badGeneIds values"<<endl;
/*
for(int i=0;i<badGeneIds.size();i++)
{
	cout<<badGeneIds[i]<<" "<<badValues[i]<<endl;
}
*/
	//real bad?
//cout<<"check is really bad"<<endl;
	for(int i=0;i<badGeneIds.size();i++)
	{
		vector<int> geneId5ps;
		vector<int> geneId3ps;

		int bdId=badGeneIds[i];
		for(int j=0;j<results.size();j++)
		{
			result_t * prt=this->getOneResult(j);
			if(prt->geneId1==bdId)
			{
				geneId3ps.push_back(prt->geneId2);
			}
			if(prt->geneId2==bdId)
			{
				geneId5ps.push_back(prt->geneId1);
			}
		}

		FocalRegionHandler frh;

		vector<region_to_map_t> vtp1;
		vector<region_to_map_t> vtup1;
		for(int x=0;x<geneId5ps.size();x++)
		{
			region_to_map_t ret;
			int gid=geneId5ps[x];
			ret.tid=g.getTid(gid);
			ret.strand=g.getStrand(gid);
			ret.lpos=g.getLimitLeft(gid);
			ret.rpos=g.getLimitRight(gid);
			vtp1.push_back(ret);
		}

		frh.getUion(vtp1,vtup1);

		vector<region_to_map_t> vtp2;
		vector<region_to_map_t> vtup2;
		for(int x=0;x<geneId3ps.size();x++)
		{
			region_to_map_t ret;
			int gid=geneId3ps[x];
			ret.tid=g.getTid(gid);
			ret.strand=g.getStrand(gid);
			ret.lpos=g.getLimitLeft(gid);
			ret.rpos=g.getLimitRight(gid);
			vtp2.push_back(ret);
		}

		frh.getUion(vtp2,vtup2);



		int realBad=vtup1.size()+vtup2.size();

		//cout<<"realbad="<<realBad<<vtup1.size()<<" "<<vtup2.size()<<endl;;

	    if(realBad>=largeNum)//para
	    {
	    	badGeneIds2.push_back(badGeneIds[i]);
	    }

	}

	for(int j=0;j<badGeneIds2.size();j++)
	{
		int bdId=badGeneIds2[j];
		for(int i=0;i<results.size();i++)
		{
			result_t * prt=this->getOneResult(i);
			if(prt->geneId1==bdId || prt->geneId2==bdId)
			{
				results[i].realPrint=0; //Do not print the result if it has a geneid in badGeneIds2
				cout<<"Removed record:"<<results[i].nm5p<<">>"<<results[i].nm3p<<endl;
//cout<<"Result "<<i<<" removed"<<endl;
			}
		}
	}
	return 0;
}
int Result::combineRecord(Gene & g)
{
//cout<<"combine"<<endl;
	if(results.size()>1)
	for(int i=0;i<results.size()-1;i++)
	{
		result_t * prtI=this->getOneResult(i);
		int en1=prtI->enrnas.size();
		int sp1=prtI->sprnas.size();
		int tidI1=g.getTid(prtI->geneId1);
		int tidI2=g.getTid(prtI->geneId2);
		int lI1=g.getLimitLeft(prtI->geneId1);
		int lI2=g.getLimitLeft(prtI->geneId2);
		for(int j=i+1;j<results.size();j++)
		{

//cout<<"ij"<<i<<" "<<j<<endl;
			result_t * prtJ=this->getOneResult(j);
			int en2=prtJ->enrnas.size();
			int sp2=prtJ->sprnas.size();
			int tidJ1=g.getTid(prtJ->geneId1);
			int tidJ2=g.getTid(prtJ->geneId2);
			if(tidI1!=tidJ1 || tidI2!=tidJ2)
				continue;
			int lJ1=g.getLimitLeft(prtJ->geneId1);
			int lJ2=g.getLimitLeft(prtJ->geneId2);
			if(en2-en1<10 && en2-en1>-10 && sp2-sp1<10 && sp2-sp1>-10)
			{
				if( lI1 < lJ1 + 1000000 && lI1 + 1000000 > lJ1 && lI2 < lJ2 + 1000000 && lI2 + 1000000 > lJ2)
				{
					int share=0;
					for(int x=0;x<prtI->sprnas.size();x++)
					{
						for(int y=0;y<prtJ->sprnas.size();y++)
						{
							if(prtI->sprnas[x].name.compare(prtJ->sprnas[y].name)==0)
							{
//cout<<"name "<<prtI->sprnas[x].name<<" "<<prtJ->sprnas[y].name<<endl;
								share++;
							}
						}
					}
					int all1=prtI->sprnas.size();
					if(share>all1*0.9)
					{
//cout<<"share"<<endl;
						if(prtI->tier<prtJ->tier) {
							prtJ->realPrint=0;
							cout<<"Combined record:"<<prtJ->nm5p<<">>"<<prtJ->nm3p<<endl;
						}
						else if(prtI->tier==prtJ->tier && prtI->realPrint!=0 && prtJ->realPrint!=0)
						{
							prtJ->realPrint=0;
							cout<<"Combined record2:"<<prtJ->nm5p<<">>"<<prtJ->nm3p<<endl;
							std::size_t found=prtI->nm5p.find(prtJ->nm5p);
//cout<<prtI->nm5p<<" "<<prtJ->nm5p<<endl;
							if(found==std::string::npos)
							{
								prtI->nm5p.append("/");
								prtI->nm5p.append(prtJ->nm5p);
//cout<<prtI->nm5p<<endl;
							}
							found=prtI->nm3p.find(prtJ->nm3p);
							if(found==std::string::npos)
							{
//cout<<prtI->nm3p<<" "<<prtJ->nm3p<<endl;
								prtI->nm3p.append("/");
								prtI->nm3p.append(prtJ->nm3p);
//cout<<prtI->nm3p<<endl;
							}


						}

					}
				}
			}
		}
	}
	return 0;
}

std::string Result::createFilename(char* dir, char* file){
	std::string filename;
	char* divider = "/";
	filename = dir;
	filename += divider;
	filename += file;
	return filename;
}


int Result::printSummary(char* filename, Gene & g, int isRunningNormal, int largeNum, char* dir) {

	sort(results.begin(),results.end(),my_sort_result_func);

	for(int i=0;i<results.size();i++)
		results[i].realPrint=1;
	

	removeMultiple(g,largeNum);
	combineRecord(g);

	std::string outfile = createFilename(dir, filename);


	ofstream outFile(outfile);
	{
		outFile<<"Fusion_Candidate\t";
		outFile<<"5_Prime\t";
		outFile<<"3_Prime\t";
		outFile<<"Reciprocal\t";
		outFile<<"FCirc\t";
		outFile<<"Tier\t";
		outFile<<"Type\t";
		outFile<<"EN_RNA\t";
		outFile<<"SP_RNA\t";
		if(indi>1)
		{

			if(isRunningNormal==0)
			{
				outFile<<"EN_DNA_T\t";
				outFile<<"SP_DNA_T\t";
			}
			else
			{
				outFile<<"EN_DNA_N\t";
				outFile<<"SP_DNA_N\t";
			}
			if(indi>2)
			{
				if(isRunningNormal==0)
				{
					outFile<<"EN_DNA_N\t";
					outFile<<"SP_DNA_N\t";
				}
				else
				{
					outFile<<"Number_1\t";
					outFile<<"Number_2\t";
					cout<<"Warning: when run with -normal, at most two data sets are needed. normal rna and normal dna."<<endl;
				}
			}
		}

        

        
        
		outFile<<"Splicings"<<endl;

		int x=0;
		for(int i=0;i<results.size();i++)
		{
			result_t rt=results[i];
			if(rt.realPrint==0)
				continue;
			outFile<<++x<<"\t"; //Candidate number
			//5' and 3' gene names
			//(Re)set in checkAllPrime()
			outFile<<rt.nm5p<<"\t"<<rt.nm3p<<"\t"; 

            ////for bk
	    //cout<<"in sum 1"<<endl;
            break_point_record_t bkt;
            bkt.nm5p=rt.nm5p;
            bkt.nm3p=rt.nm3p;
            //bkvec.push_back(bkt);
            ////
            //cout<<"get"<<bkt.nm5p<<" "<<bkt.nm3p<<"and pushed"<<endl;
            
			if(rt.isReci==1)
				outFile<<"Y\t";
			else
				outFile<<"N\t";
			if(rt.hasFcirc==1)
				outFile<<"Y\t";
			else
				outFile<<"N\t";


			outFile<<rt.tier<<"\t";

	    //// for bk 2
            bkt.tier=rt.tier;
	
			int t1=6;

			for(int x=0;x<rt.types.size();x++)
			{
				if(rt.types[x]<t1)
				{
					t1=rt.types[x];
				}
			}
			
			outFile<<getType(t1)<<"\t"; // Type of fusion
             //// for bk 2
             if(t1==2)
				bkt.isRT==1;
	     	else
				bkt.isRT==0;
             bkvec.push_back(bkt);  
			//Read counts
			outFile<<rt.numOfEnRna<<"\t"<<rt.numOfSpRna<<"\t";
			if(indi > 1)
			{
				outFile<<rt.numOfEnDnaT<<"\t"<<rt.numOfSpDnaT<<"\t";
				if(indi>2)
				{
					outFile<<rt.numOfEnDnaN<<"\t"<<rt.numOfSpDnaN<<"\t";
				}
			}
			for(int x=0;x<rt.types.size();x++)
			{
				if(rt.sameDirAsMain[x]==1){
					outFile<<rt.nm5p<<">>"<<rt.nm3p;
				} else {
					outFile<<rt.nm3p<<">>"<<rt.nm5p;
				}
				/*
				if(rt.primeOKs[x]==1)
				{
					outFile<<rt.nm5p<<">>"<<rt.nm3p;
				}
				else
				{
					outFile<<rt.nm3p<<">>"<<rt.nm5p;
				}
				*/
				outFile<<"("<<getType(rt.types[x])<<" ";
				if(rt.canos[x]==1)
					outFile<<"Canonical ";
				outFile<<rt.numOfsps[x]<<");";
			}

			outFile<<endl;

		}
	}

	return 0;
}

vector<int> Result::getAllSprnas(result_t * prt, vector<split_rna_t> & sts)
{

	vector<int> support;
	int sum=0;
	//iterate over junctions
	for(int i=0;i<prt->numOfsps.size();i++) 
	{
		//if(prt->canos[i]==1) ??
		sts.push_back(prt->sprnas[sum]); // Collect a sprna for the junction
		support.push_back(prt->numOfsps[i]); // Collect the # of sprna at the junction
		sum+=prt->numOfsps[i]; // Increase sum for next round of collection
	}

	return support;
}



void Result::checkMultiJunction(int index, Gene& g) 
{

	result_t * prt = this->getOneResult(index); //Grab fusion

	//Set defaults
	prt->isReci = 0;
	prt->hasFcirc = 0;

	vector<split_rna_t> sts;
	vector<int> support = getAllSprnas(prt, sts); // Grab all split RNA reads
	//support == # supporting reads for each junction
	//sts describes each junction being supported

	// identify the junction with the most support
	// Assume it is the real fusion
	int max_support = 0;
	int index_of_max = 0;
	for(int i=0;i<support.size();i++)
	{
		if(support[i]>max_support)
		{
			max_support = support[i];
			index_of_max = i;
		}
	}

	split_rna_t main_junction = sts[index_of_max];
	//// identify the exons associated with the primary junction
	//First exon
	int p1;
	if(main_junction.bkLeft1==1)
		p1 = main_junction.pos1;
	else
		p1 = main_junction.pos1+main_junction.len1-1;
	

	int is5p, tid, strand, pos1, pos2, main_exonNum;
	string main_name;
	int tranId;
	
	g.getBestExon(main_junction.geneId1, p1, 
		main_junction.bkLeft1, is5p, tid, strand,
		pos1, pos2, main_name, main_exonNum, tranId);
	
	//Second exon
	int p2;
	if(main_junction.bkLeft2==1)
		p2 = main_junction.pos2;
	else	 
		p2 =main_junction.pos2+main_junction.len2-1;
	int is5p_2, tid_2, strand_2, pos1_2, pos2_2, main_exonNum_2;
	string main_name_2;
	int tranId_2;

	g.getBestExon(main_junction.geneId2, p2, main_junction.bkLeft2,
		is5p_2, tid_2, strand_2, 
		pos1_2, pos2_2, main_name_2, main_exonNum_2, tranId_2);
	
	string main_gname = g.getName2(main_junction.geneId1); //Was st.
	string main_gname_2 = g.getName2(main_junction.geneId2);

	
	//Label 5p and 3p gene for orientation comparisons
	string gene5p, gene3p;
	int exon5p, exon3p;
	int pos5p, pos3p;
	int strand5p, strand3p;
	int tranId5p, tranId3p;
	string tranName5p, tranName3p;
	if(is5p==1) {
		gene5p = main_gname;
		exon5p = main_exonNum;
		gene3p = main_gname_2;
		exon3p = main_exonNum_2;
		pos5p = p1;
		pos3p = p2;
		strand5p = strand;
		strand3p = strand_2;
		tranId5p = tranId;
		tranId3p = tranId_2;
		tranName5p = main_name;
		tranName3p = main_name_2;
		if (prt->nm5p!=gene5p) {
			prt->nm5p = gene5p;
			prt->nm3p = gene3p;
		}
	} else {
		gene5p = main_gname_2;
		exon5p = main_exonNum_2;
		gene3p = main_gname;
		exon3p = main_exonNum;
		pos5p = p2;
		pos3p = p1;
		strand5p = strand_2;
		strand3p = strand;
		tranId5p = tranId_2;
		tranId3p = tranId;
		tranName5p = main_name_2;
		tranName3p = main_name;
		if (prt->nm5p!=gene5p) {
			prt->nm5p = gene5p;
			prt->nm3p = gene3p;
		}
	}
	


	//Compare all other junctions to the primary junction

	//To help keep things straight, remove magic numbers
	int inter_chrom = 0;
	int intra_chrom = 1;
	int readt = 2;
	int splice_v = 3;
	int backsplice = 4;
	int other = 5;
	vector<int> indicesOfPotentialFCircOfReci;
	vector<int> indicesOfPotentialAltsOfReci;
	int maxSupportForAlts = 0;
	
	bool reciWasFound = false;

	for(int i=0;i<sts.size();i++)
	{
		if(i==index_of_max) {
			prt->sameDirAsMain.push_back(1);
		}
		else 
		{

			split_rna_t st = sts[i];

			string gname = g.getName2(st.geneId1);
			string gname_2 = g.getName2(st.geneId2);

			//Get best exons for given junction
			//First exon
			int p1;
			if(st.bkLeft1==1)
				p1=st.pos1;
			else
				p1=st.pos1+st.len1-1;

			//Second exon
			int p2;
			if(st.bkLeft2==1)
				p2=st.pos2;
			else
				p2=st.pos2+st.len2-1;

			if(!g.posIsValid(st.geneId1, p1)) {
					int tmp = p1;
					p1 = p2;
					p2 = tmp;
				}	


			//Just reuse old variables
			//int is5p, tid, strand, pos1, pos2;
			string name1;
			int exonNum;
			g.getBestExon3(st.geneId1, p1, st.bkLeft1,
				is5p, tid, strand, pos1, pos2, name1, exonNum,
				 gene5p, gname, tranId5p, tranId3p, tranId);
			
			

			//Just reuse old variables
			//int is5p_2, tid_2, strand_2, pos1_2, pos2_2;
			string name_2;
			int exonNum_2;
			g.getBestExon3(st.geneId2, p2, st.bkLeft2,
				is5p_2, tid_2, strand_2, pos1_2, pos2_2, name_2, exonNum_2,
				 gene5p, gname_2, tranId5p, tranId3p, tranId_2);

	

			string currGene5p, currGene3p;
			int currExon5p, currExon3p;
			int currPos5p, currPos3p;
			int currStrand5p, currStrand3p;
			int currTranId5p, currTranId3p;
			string currTranName5p, currTranName3p;
			//Flip as needed
			if(is5p==1) {
				currGene5p = gname;
				currExon5p = exonNum;
				currGene3p = gname_2;
				currExon3p = exonNum_2;
				currPos5p = p1;
				currPos3p = p2;
				currStrand5p = strand;
				currStrand3p = strand_2;
				currTranId5p = tranId;
				currTranId3p = tranId_2;
				currTranName5p = name1;
				currTranName3p = name_2;
				//Double check that pos is correct
				//st.geneId1->currGene5p
				//Make sure currPos5p is from that gene
				//Do not have a good way to confirm exonNum
				if(!g.posIsValid(st.geneId1, currPos5p)) {
					int tmp = currPos5p;
					currPos5p = currPos3p;
					currPos3p = tmp;
				}	
			} else {
				currGene5p = gname_2;
				currExon5p = exonNum_2;
				currGene3p = gname;
				currExon3p = exonNum;
				currPos5p = p2;
				currPos3p = p1;
				currStrand5p = strand_2;
				currStrand3p = strand;
				currTranId5p = tranId_2;
				currTranId3p = tranId;
				currTranName5p = name_2;
				currTranName3p = name1;
				//st.geneId1->gname->currGene3p
				//Make sure currPos3p is from that gene
				//Do not have a good way to confirm exonNum
				if(!g.posIsValid(st.geneId1, currPos3p)) {
					int tmp = currPos5p;
					currPos5p = currPos3p;
					currPos3p = tmp;
				}
			}

			


			
			//Set info for printing later
			if(currGene5p==gene5p) {
				prt->sameDirAsMain.push_back(1);
			} else {
				prt->sameDirAsMain.push_back(0);
			}


			//When comparing exon numbers, ensure they are from the correct transcripts,
			//else use pos information.
			//Relying on exon number helps remove junctions that are just a couple bases
			//different from each other, so don't want to rely solely on pos if possible
			if(currGene5p==gene5p)
			{
				if(currExon5p<=exon5p && currExon3p>=exon3p && currTranName5p==tranName5p && currTranName3p==tranName3p) {
					//Almost certainly a splice variant junction
					//unless the main junction is an alt splice that is more common
					//than the real fusion
					results[index].types[i] = splice_v; 
					//Shorter, linear transcript than main fusion
				}
				else if(currExon5p>exon5p && currExon3p<exon3p && currTranName5p==tranName5p && currTranName3p==tranName3p) {
					//Either the main function is a highly common alt splice
					//of this real fusion, or it is an fcirc of a reciprocal fusion (unlikely)
					indicesOfPotentialFCircOfReci.push_back(i); //Add to list to check later
					//prt->isReci = 1;
					//prt->hasFcirc = 1;
					//results[index].types[i] = 4; 
					//Other option is that the main junction was an alternative splice,
					//with more read support than the actual fusion
				}
				else { //Attempt to annotate based on pos as last resort
					//currGene5p == gene5p
					//strand5p == strand of main 5' gene
					//strand3p == strand of currStrand3p == strand of main 3' gene
					if(strand5p==0 && strand3p==0) {
						if(currPos5p > pos5p && currPos3p < pos3p) {
							//Also possible that this is the real fusion,
							//and the main fusion is an alt splice of it
							indicesOfPotentialFCircOfReci.push_back(i);
							//prt->isReci = 1;
							//prt->hasFcirc = 1;
							//results[index].types[i] = 4; 
						} else if (currPos5p < pos5p && currPos3p > pos3p) {
							results[index].types[i] = splice_v;
						} else {
							results[index].types[i] = other;
						}
					} else if (strand5p==0 && strand3p==1) {
						if(currPos5p > pos5p && currPos3p > pos3p) {
							indicesOfPotentialFCircOfReci.push_back(i);
							//prt->isReci = 1;
							//prt->hasFcirc = 1;
							//results[index].types[i] = 4;
						} else if (currPos5p < pos5p && currPos3p < pos3p) {
							results[index].types[i] = splice_v;
						} else {
							results[index].types[i] = other;
						}
					} else if (strand5p==1 && strand3p==0) {
						if(currPos5p < pos5p && currPos3p < pos3p) {
							indicesOfPotentialFCircOfReci.push_back(i);
							//prt->isReci = 1;
							//prt->hasFcirc = 1;
							//results[index].types[i] = 4;
						} else if (currPos5p > pos5p && currPos3p > pos3p) {
							results[index].types[i] = splice_v;
						} else {
							results[index].types[i] = other;
						}
					} else if (strand5p==1 && strand3p==1) {
						if(currPos5p < pos5p && currPos3p > pos3p) {
							indicesOfPotentialFCircOfReci.push_back(i);
							//prt->isReci = 1;
							//prt->hasFcirc = 1;
							//results[index].types[i] = 4;
						} else if (currPos5p > pos5p && currPos3p < pos3p) {
							results[index].types[i] = splice_v;
						} else {
							results[index].types[i] = other;
						}
					}
				}
			}
			else if(currGene5p==gene3p)
			{
				if(currExon5p==(exon3p-1) && currExon3p==(exon5p+1) && 
					currTranName5p==tranName3p && currTranName3p==tranName5p)
				{
					prt->isReci = 1;
					reciWasFound = true;
				}
				else if(currExon5p>exon3p && currExon3p<exon5p && 
					currTranName5p==tranName3p && currTranName3p==tranName5p)
				{
					prt->hasFcirc = 1;
					results[index].types[i] = backsplice;
				}
				else if(currExon5p<=(exon3p-1) && currExon3p>=(exon5p+1) &&
					currTranName5p==tranName3p && currTranName3p==tranName5p)
				{
					indicesOfPotentialAltsOfReci.push_back(i);
					if(support[i] > maxSupportForAlts)
						maxSupportForAlts = support[i];
					//prt->isReci = 1;
					//results[index].types[i] = 3;
				}
				else {
					//currGene5p == gene3p
					//strand5p = strand of main 5' gene
					//currStrand5p = strand of main 3' gene
					//Unlikely to be able to call reciprocal fusion based on pos
					if(strand5p==0 && currStrand5p==0) {
						if(currPos5p > pos3p && currPos3p < pos5p) { //FCirc
							prt->hasFcirc = 1;
							results[index].types[i] = backsplice;
						} else if (currPos5p < pos3p && currPos3p > pos5p) { //Alt of reci
							//Potentially the actual reciprocal fusion
							indicesOfPotentialAltsOfReci.push_back(i);
							if(support[i] > maxSupportForAlts)
								maxSupportForAlts = support[i];
							//prt->isReci = 1;
							//results[index].types[i] = 3;
						} else {
							results[index].types[i] = other;
						}
					} else if (strand5p==0 && currStrand5p==1) {
						if(currPos5p < pos3p && currPos3p < pos5p) {
							prt->hasFcirc = 1;
							results[index].types[i] = backsplice;
						} else if (currPos5p > pos3p && currPos3p > pos5p) { //Alt of reci
							indicesOfPotentialAltsOfReci.push_back(i);
							if(support[i] > maxSupportForAlts)
								maxSupportForAlts = support[i];
							//prt->isReci = 1;
							//results[index].types[i] = 3;
						} else {
							results[index].types[i] = other;
						}
					} else if (strand5p==1 && currStrand5p==0) {
						if(currPos5p > pos3p && currPos3p > pos5p) {
							prt->hasFcirc = 1;
							results[index].types[i] = backsplice;
						} else if (currPos5p < pos3p && currPos3p < pos5p) { //Alt of reci
							indicesOfPotentialAltsOfReci.push_back(i);
							if(support[i] > maxSupportForAlts)
								maxSupportForAlts = support[i];
							//prt->isReci = 1;
							//results[index].types[i] = 3;
						} else {
							results[index].types[i] = other;
						}
					} else if (strand5p==1 && currStrand5p==1) {
						if(currPos5p < pos3p && currPos3p > pos5p) {
							prt->hasFcirc = 1;
							results[index].types[i] = backsplice;
						} else if(currPos5p > pos3p && currPos3p < pos5p) { //Alt of reci
							indicesOfPotentialAltsOfReci.push_back(i);
							if(support[i] > maxSupportForAlts)
								maxSupportForAlts = support[i];
							//prt->isReci = 1;
							//results[index].types[i] = 3;
						} else {
							results[index].types[i] = other;
						}
					}
				}
			}

		}
	}

	//Revisit this once support for DNA has been added. 
	//Expect DNA support for fusions but not alt splices

	//If a reciprocal was found, these values are likely alternative splices
	//If a reciprocal was not found (based on exons), these may be the real reciprocal (based on pos)
	for(int x=0;x<indicesOfPotentialAltsOfReci.size();x++){
		if(reciWasFound){
			results[index].types[indicesOfPotentialAltsOfReci[x]] = splice_v;
			prt->isReci = 1; //Should already be 1
		} else {
			//Reciprocal was not found, but these look like reciprocals
			//Leave types info the same as before (Inter or Intra chromosomal fusion)
			//but update isReci
			if(support[x]!=maxSupportForAlts) {
				//Multiple potential reci were found. This has less support,
				//so treat it as the alt. Leave the type of the max alone
				results[index].types[indicesOfPotentialAltsOfReci[x]] = splice_v;
			}
			prt->isReci = 1;
		}
	}

	for(int x=0;x<indicesOfPotentialFCircOfReci.size();x++) {
		if(reciWasFound){
			//If a reciprocal was found (based on exon counts)
			//then treat these as fcircRNAs of the reciprocal fusion
			results[index].types[indicesOfPotentialFCircOfReci[x]] = backsplice;
		} else {
			//If no reciprocal was found (based on exon counts)
			//then assume this is just an (odd) alternative splice of the main fusion
			//(potentially even the main fusion)
			results[index].types[indicesOfPotentialFCircOfReci[x]] = splice_v;
		}
	}




	return;

}

int Result::checkALLPrime(Gene & g) {

	for(int i=0;i<results.size();i++)
	{
		int ok=0;
		int all=results[i].primeOKs.size();
		for(int j=0;j<all;j++)
		{
			if(results[i].primeOKs[j]==1)
				ok++;
		}

		// ok==all if all primeOKs==1. This would mean 5' and 3' gene are labelled correctly
		//if ok<all*0.5, then the gene labeled as 3' is more likely to be the 5' gene (got it reversed),
		//because >50% of junctions are in the other direction.
		//ok should either ==all or ==0 unless a reciprocal fusion or fcircRNA is present
		
		if(ok<all*0.5)
		{

			int tmp=results[i].geneId1;
			results[i].geneId1=results[i].geneId2;
			results[i].geneId2=tmp;

			string name=results[i].nm5p;
			results[i].nm5p=results[i].nm3p;
			results[i].nm3p=name;

			

			for(int j=0;j<all;j++)
			{
				results[i].primeOKs[j]=1-results[i].primeOKs[j];
			}
		}


		// Set isReci and isFcirc values
		// Assumes that the most well supported junction is the real fusion
		// We don't really know enough about fcircRNA to know if this is a valid assumption
		

		checkMultiJunction(i,g);
		
		/*
		//Method for setting isReci from original INTEGRATE tool
		if(ok!=all && ok!=0)
		{
			results[i].isReci=1;
		}
		else
		{
			results[i].isReci=0;
		}
		*/


	}


	for(int i=0;i<results.size();i++)
	{
		int cano=0;
		int all=results[i].primeOKs.size();
		for(int j=0;j<all;j++)
		{
			if(results[i].canos[j]==1)
			{
				cano=1;
				break;
			}

		}
		if(cano==1)
		{
			results[i].isCanonical=1;
		}
		else
		{
			results[i].isCanonical=0;
		}

	}


	return 0;
}



int Result::getSize() {
	return results.size();
}

//Dec 7,2015, add for printAllJunctions
int bestSprnaVec(result_t * prt,vector<split_rna_t> & st_vc)
{
    int sum=0;
    for(int i=0;i<prt->numOfsps.size();i++)
    {
        if(prt->canos[i]==1)
        {
            split_rna_t st=prt->sprnas[sum];
            st_vc.push_back(st);
        }
        sum+=prt->numOfsps[i];
    }
    return 0;
}


int bestSprna(result_t * prt,split_rna_t & st)
{
	int num=0;
	int sum=0;
	for(int i=0;i<prt->numOfsps.size();i++)
	{
		if(prt->numOfsps[i]>num && prt->canos[i]==1)
		{
			num=prt->numOfsps[i];
			st=prt->sprnas[sum];
		}
		sum+=prt->numOfsps[i];
	}
	return 0;
}

int bestSprnaGen(result_t * prt,split_rna_t & st)
{
	int num=0;
	int sum=0;
	for(int i=0;i<prt->numOfsps.size();i++)
	{
		if(prt->numOfsps[i]>num)
		{
			num=prt->numOfsps[i];
			st=prt->sprnas[sum];
		}
		sum+=prt->numOfsps[i];
	}
	return 0;
}


int bestSpdna(result_t * prt,split_dna_t & st)
{
    
    int best=0;
    int maxmin=0;
    for(int i=0;i<prt->spdna1.size();i++)
    {
        int len1=prt->spdna1[i].len1;
        int len2=prt->spdna1[i].len2;
        int minlen=len1;
        if(len2<len1)
            minlen=len2;

        if (minlen>maxmin) {
            maxmin=minlen;
            best=i;
        }
    }
    st=prt->spdna1[best];
	return 0;
}


int RefPrinter(Reference &ref, int tid, uint32_t aa, uint32_t bb, int strand, int is5p, ofstream & outFile)
{


	        uint32_t refaa=ref.to_ref_pos(tid,aa);
	        uint32_t refbb=ref.to_ref_pos(tid,bb);
	        for(uint32_t x=refaa;x<=refbb;x++)
	        	outFile<<ref.getRefChar(x);
	        outFile<<"\t";

	        if(refbb-refaa+1>=150)
	        {
	                if((strand==0 && is5p==1) || (strand==1 && is5p==0))
	                {
	                        for(uint32_t x=refbb-150+1;x<=refbb;x++)
	                        	outFile<<ref.getRefChar(x);
	                }
	                else
	                {
	                        for(uint32_t x=refaa;x<=refaa+150-1;x++)
	                        	outFile<<ref.getRefChar(x);
	                }
	        }
	        else
	        	outFile<<"NA";

	        return 0;
}





//bkfileBEDPE is deprecated
int Result::printExons(char* file, Gene& g, Reference & ref, int isRunningNormal, char * bkfile, char * bkfileBEDPE, char * bkfileVCF, char * refname, char * sample_name, char* dir) {

	std::string filename = createFilename(dir,file);
	ofstream outFile(filename);

	outFile<<"#Id\t5p\t3P\t5P_Transcipt\t5P_Exon\t5P_Strand\t5P_Exon_Chr\t5P_Exon_Start\t5P_Exon_End\t5P_Exon_Seq\t5P_Exon_150\t3P_Transcript\t3P_Exon\t3P_Exon_Strand\t3P_Exon_Chr\t3P_Exon_Start\t3P_Exon_END\t3P_Exon_Seq\t3P_Exon_150"<<endl;

    int index=0;
    
	for(int i=0;i<results.size();i++)
	{
		result_t * prt=this->getOneResult(i);
		if(prt->tier>3)
			break;
		if(prt->realPrint==0)
            continue;

        
        

		outFile<<index+1<<"\t"; //ID

		split_rna_t st;
		bestSprna(prt,st);

	    int p1;
	    if(st.bkLeft1==1)
	        p1=st.pos1;
	    else
	        p1=st.pos1+st.len1-1;


	    int p2;
	    if(st.bkLeft2==1)
	        p2=st.pos2;
	    else
	        p2=st.pos2+st.len2-1;

		int is5p;
		int tid;
		int strand;
		int pos1;
		int pos2;
		string name;
		int exonNum;
		int tranId;

		g.getBestExon(st.geneId1,p1, st.bkLeft1,
				is5p, tid, strand, pos1, pos2, name, exonNum, tranId);

		int is5p_2;
		int tid_2;
		int strand_2;
		int pos1_2;
		int pos2_2;
		string name_2;
		int exonNum_2;

		g.getBestExon(st.geneId2,p2, st.bkLeft2,
				is5p_2, tid_2, strand_2, pos1_2, pos2_2, name_2, exonNum_2, tranId);



        
        
		if(is5p==1 && is5p_2==0)
		{

			outFile<<g.getName2(st.geneId1)<<"\t"<<g.getName2(st.geneId2)<<"\t";
			outFile<<name<<"\t"<<exonNum<<"\t";
			if(strand==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(tid)<<"\t";
			outFile<<pos1<<"\t"<<pos2<<"\t";
			RefPrinter(ref,tid,pos1,pos2,strand,is5p,outFile);
			outFile<<"\t";
			outFile<<name_2<<"\t"<<exonNum_2<<"\t";
			if(strand_2==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(tid_2)<<"\t";
			outFile<<pos1_2<<"\t"<<pos2_2<<"\t";
			RefPrinter(ref,tid_2,pos1_2,pos2_2,strand_2,is5p_2,outFile);
			outFile<<endl;
            
            ////for bk
            //cout<<"in exon 1"<<endl;
            bkvec[index].tid1=tid;
            bkvec[index].tid2=tid_2;
            bkvec[index].isExon=1;
            bkvec[index].swp=0;            

            if(strand==0)
            {
                bkvec[index].exonbk1=pos2;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand1=0;
            }
            else
            {
                bkvec[index].exonbk1=pos1;
                bkvec[index].seqLeft1=0;
		bkvec[index].gStrand1=1;
            }
            
            if(strand_2==0)
            {
                bkvec[index].exonbk2=pos1_2;
                bkvec[index].seqLeft2=0;
		bkvec[index].gStrand2=0;
            }
            else
            {
                bkvec[index].exonbk2=pos2_2;
                bkvec[index].seqLeft2=1;
		bkvec[index].gStrand2=1;
            }
            bkvec[index].splitrna=st;
            //////
		//cout<<"left 1"<<endl;
            
		}
		else
		{
			outFile<<g.getName2(st.geneId2)<<"\t"<<g.getName2(st.geneId1)<<"\t";
			outFile<<name_2<<"\t"<<exonNum_2<<"\t";
			if(strand_2==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(tid_2)<<"\t";
			outFile<<pos1_2<<"\t"<<pos2_2<<"\t";
			RefPrinter(ref,tid_2,pos1_2,pos2_2,strand_2,is5p_2,outFile);
			outFile<<"\t";
			outFile<<name<<"\t"<<exonNum<<"\t";
			if(strand==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(tid)<<"\t";
			outFile<<pos1<<"\t"<<pos2<<"\t";
			RefPrinter(ref,tid,pos1,pos2,strand,is5p,outFile);
			outFile<<endl;
            
            ////for bk
            //cout<<"in exon 2"<<endl; 
            bkvec[index].tid1=tid_2;
            bkvec[index].tid2=tid;
            bkvec[index].isExon=1;
            bkvec[index].swp=1;

            if(strand_2==0)
            {
                bkvec[index].exonbk1=pos2_2;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand1=0;
            }
            else
            {
                bkvec[index].exonbk1=pos1_2;
                bkvec[index].seqLeft1=0;
		bkvec[index].gStrand1=1;
            }
            if(strand==0)
            {
                bkvec[index].exonbk2=pos1;
                bkvec[index].seqLeft2=0;
		bkvec[index].gStrand2=0;
            }else
            {
                bkvec[index].exonbk2=pos2;
                bkvec[index].seqLeft2=1;
		bkvec[index].gStrand2=1;
            }
            bkvec[index].splitrna=st;
            /////
	//cout<<"%%%%%%"<<bkvec[index].seqLeft1<<" "<<bkvec[index].seqLeft2<<endl;
	    //cout<<"left 2"<<endl;
		}

        if(prt->tier==1)
        {
	////
	    //cout<<"get best dna"<<endl;
            bkvec[index].rna_only=0;
            split_dna_t sdt;
            bestSpdna(prt,sdt);
            bkvec[index].splitdna=sdt;
	    //cout<<"left here"<<endl;
        }
        else
        {
            bkvec[index].rna_only=1;
        }
        
        index++;
        
	}
    
    /////// for bk non canonical
    //index=0;
    for(int i=0;i<results.size();i++)
	{
		result_t * prt=this->getOneResult(i);
		if(prt->tier<=3)
			continue;
		if(prt->realPrint==0)
            continue;
        
        
		split_rna_t st;
		bestSprnaGen(prt,st);
        
	    int p1;
	    if(st.bkLeft1==1)
	        p1=st.pos1;
	    else
	        p1=st.pos1+st.len1-1;
        
        
	    int p2;
	    if(st.bkLeft2==1)
	        p2=st.pos2;
	    else
	        p2=st.pos2+st.len2-1;

        int is5p;
		int tid;
		int strand;
        g.getStrandnPrimenTid(st.geneId1,st.bkLeft1, is5p,tid,strand);

		int is5p_2;
		int tid_2;
		int strand_2;
        
        g.getStrandnPrimenTid(st.geneId2,st.bkLeft2, is5p_2,tid_2, strand_2);

        
		if(is5p==1 && is5p_2==0)
		{
           
           ////for bk

		//cout<<"in non cano 1"<<endl; 
            bkvec[index].tid1=tid;
            bkvec[index].tid2=tid_2;
            bkvec[index].isExon=0;
            bkvec[index].swp=0;

		//cout<<"tid,tid2,isExon"<<tid<<","<<tid_2<<","<<0<<endl;		

		//cout<<"strand="<<strand<<endl;

            if(strand==0)
            {
		//cout<<"if1"<<endl;
                bkvec[index].exonbk1=0;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand1=0;
            }
            else
            {
		//cout<<"else1"<<endl;
                bkvec[index].exonbk1=0;
                bkvec[index].seqLeft1=0;
		bkvec[index].gStrand1=1;
            }
            //cout<<"strand_2="<<strand_2<<endl;
            if(strand_2==0)
            {
		//cout<<"if2"<<endl;
                bkvec[index].exonbk2=0;
                bkvec[index].seqLeft2=0;
		bkvec[index].gStrand2=0;
            }
            else
            {
		//cout<<"else2"<<endl;
                bkvec[index].exonbk2=0;
                bkvec[index].seqLeft2=1;
		bkvec[index].gStrand2=1;
            }
		//cout<<"prepare to copy"<<endl;
            bkvec[index].splitrna=st;
            //cout<<"left non cano 1"<<endl;
            
		}
		else
		{

            
            ////for bk
           
		//cout<<"in exon non-cano 2"<<endl; 
            bkvec[index].tid1=tid_2;
            bkvec[index].tid2=tid;
            bkvec[index].isExon=0;
            bkvec[index].swp=1;           


             

 
            if(strand_2==0)
            {
                bkvec[index].exonbk1=0;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand1=0;
            }
            else
            {
                bkvec[index].exonbk1=0;
                bkvec[index].seqLeft1=0;
		bkvec[index].gStrand1=1;
            }
            if(strand==0)
            {
                bkvec[index].exonbk2=0;
                bkvec[index].seqLeft2=0;
		bkvec[index].gStrand2=0;
            }else
            {
                bkvec[index].exonbk2=0;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand2=1;
            }
            bkvec[index].splitrna=st;
            
		//cout<<"left non cano 2"<<endl;

		}
        if(prt->tier==4)
        {
	//cout<<"non cano dna"<<endl;
            bkvec[index].rna_only=0;
            split_dna_t sdt;
            bestSpdna(prt,sdt);
            bkvec[index].splitdna=sdt;
	//cout<<"left here"<<endl;
        }
        else
        {
            bkvec[index].rna_only=1;
        }
        index++;
	}
    
	//cout<<"call bk"<<endl;
    	BreakPoint bkobj;
    	bkobj.getBreakPoints(bkvec,bkfile,bkfileBEDPE,bkfileVCF,refname,ref,sample_name,dir);
    	//cout<<"after bk"<<endl;
    
	return 0;
}




int getRefExonSeq(Reference &ref, int tid, uint32_t aa, uint32_t bb, int strand, vector<char> & seq)
{
    
    uint32_t refaa=ref.to_ref_pos(tid,aa);
    uint32_t refbb=ref.to_ref_pos(tid,bb);
    if(strand==0)
    {
        for(uint32_t x=refaa;x<=refbb;x++)
            seq.push_back(ref.getRefChar(x));
    }
    else
    {
        for(uint32_t x=refbb;x>=refaa;x--)
            seq.push_back(getCharComp(ref.getRefChar(x)));
    }
    
    return 0;
}


typedef struct
{
    int fusion_id;
    junction_t p5;
    junction_t p3;
    vector<char> seq1;
    vector<char> seq2;
    int isInframe;
	int type;
} fusion_junction_t;

typedef struct
{
	int curr_id;
	string name;
	int fusion5strand;
	int fusion3strand;
	uint32_t fusion5breakPos;
	uint32_t fusion3breakPos;
	string fusion5chrm;
	string fusion3chrm;
	uint32_t fusion5backPos;
	uint32_t fusion3backPos;
	string gene5;
	string gene3;
	int hasFcirc;
	//int mainNumSpReads;
	//int backspliceSpReads;
} fcirc_t;



// have repeats
vector<fusion_junction_t> fjtvec;

int Result::getAllJunctionsStep1(Gene& g, Reference & ref) {
        
    return 0;
}


int getInframe(fusion_junction_t & ft)
{
//cout<<"######"<<ft.p5.isCoding<<" "<<ft.p3.isCoding<<" "<<ft.p5.coding_left<<" "<<ft.p3.coding_left<<endl;
    if (ft.p5.isCoding==1 && ft.p3.isCoding==1 && ft.p5.coding_left>=0 && ft.p3.coding_left>=0)
    {
        if((ft.p5.coding_left+ft.p3.coding_left)%3==0)
        {
            return 2;
        }
    }
    
    if((ft.p5.isCoding==0 || ft.p5.coding_left<0 ) && ft.p3.isCoding==1 && (ft.p3.coding_left==-1 || ft.p3.coding_left==0))//this has to be ==-1
    {
        return 1;
    }
    return 0;
}

bool my_sort_for_in_frame(fusion_junction_t i, fusion_junction_t j) // id, junc5p, junc3p,in-frame, short, short
{
    if(i.fusion_id<j.fusion_id)
    {
        return true;
    }
    else if(i.fusion_id==j.fusion_id)
    {
        int p5_pos1;
        if(i.p5.strand==0)
            p5_pos1=i.p5.pos2;
        else
            p5_pos1=i.p5.pos1;
                
        int p5_pos2;
        if(j.p5.strand==0)
            p5_pos2=j.p5.pos2;
        else
            p5_pos2=j.p5.pos1;
                
                
                if(p5_pos1<p5_pos2)
                {
                    return true;
                }
                else if(p5_pos1==p5_pos2)
                {
                    
                    int p3_pos1;
                    if(i.p3.strand==0)
                        p3_pos1=i.p3.pos1;
                    else
                        p3_pos1=i.p3.pos2;
                            
                    int p3_pos2;
                    if(j.p3.strand==0)
                        p3_pos2=j.p3.pos1;
                    else
                        p3_pos2=j.p3.pos2;
                            
                    if(p3_pos1<p3_pos2)
                    {
                        return true;
                    }
                    else if(p3_pos1==p3_pos2)
                    {
                        if(i.isInframe>j.isInframe)
                        {
                            return true;
                        }
                        else if(i.isInframe==j.isInframe)
                        {
                            if(i.p5.pos2-i.p5.pos1 < j.p5.pos2-j.p5.pos1)
                            {
                                return true;
                            }
                            else if(i.p5.pos2-i.p5.pos1 == j.p5.pos2-j.p5.pos1)
                            {
                                if(i.p3.pos2-i.p3.pos1 < j.p3.pos2-j.p3.pos1)
                                {
                                    return true;
                                }
                                else
                                {
                                    return false;
                                }
                            }
                            else
                            {
                                return false;
                            }
                        }
                        else
                        {
                            return false;
                        }
                    }
                    else
                    {
                        return false;
                    }
                    
                }
        else
        {
            return false;
        }
        
    }
    else
    {
        return false;
    }
    
}



int Result::getAllJunctionsStep2(char * filename, Gene& g, Reference & ref) {

    
    return 0;


}

bool my_sort_for_peptide(fusion_junction_t i, fusion_junction_t j) // id, junc5p, junc3p,!=-1, !=-1, base, base, short, short
{
    if(i.fusion_id<j.fusion_id)
    {
        return true;
    }
    else if(i.fusion_id==j.fusion_id)
    {
        int p5_pos1;
        if(i.p5.strand==0)
            p5_pos1=i.p5.pos2;
        else
            p5_pos1=i.p5.pos1;
        
        int p5_pos2;
        if(j.p5.strand==0)
            p5_pos2=j.p5.pos2;
        else
            p5_pos2=j.p5.pos1;
        
        
        if(p5_pos1<p5_pos2)
        {
            return true;
        }
        else if(p5_pos1==p5_pos2)
        {
            
            int p3_pos1;
            if(i.p3.strand==0)
                p3_pos1=i.p3.pos1;
            else
                p3_pos1=i.p3.pos2;
            
            int p3_pos2;
            if(j.p3.strand==0)
                p3_pos2=j.p3.pos1;
            else
                p3_pos2=j.p3.pos2;
            
            if(p3_pos1<p3_pos2)
            {
                return true;
            }
            else if(p3_pos1==p3_pos2)
            {
                if(i.p5.coding_left > j.p5.coding_left)
                {
                    return true;
                }
                else if(i.p5.coding_left==j.p5.coding_left)
                {
                    if(i.p3.coding_left>j.p3.coding_left)
                    {
                        return true;
                    }
                    else if(i.p3.coding_left==j.p3.coding_left)
                    {
                        if(i.p5.pos2-i.p5.pos1 < j.p5.pos2-j.p5.pos1)
                        {
                            return true;
                        }
                        else if(i.p5.pos2-i.p5.pos1 == j.p5.pos2-j.p5.pos1)
                        {
                            if(i.p3.pos2-i.p3.pos1 < j.p3.pos2-j.p3.pos1)
                            {
                                return true;
                            }
                            else
                            {
                                return false;
                            }
                        }
                        else
                        {
                            return false;
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
            
        }
        else
        {
            return false;
        }
        
    }
    else
    {
        return false;
    }
    
}







int Result::getAllJunctionsStep3(char * filename, Gene& g, Reference & ref)
{
    
    return 0;
    
}





//SMC-RNA

typedef struct
{
    int fusion_id;
    junction_t p5;
    junction_t p3;
    vector<char> seq1;
    vector<char> seq2;
    int isInframe;
    
    int tier;
    int isCano;
    int numSpReads;
	int numEnReads;
    
    int geneId5p;
    int geneId3p;
    uint32_t pos5p;
    uint32_t pos3p;

    split_rna_t st; 

	int type; //New

} fusion_junction_2_t;




vector<fusion_junction_2_t> fjtvec2;


int bestSprnaVec2(result_t * prt,vector<split_rna_t> & st_vc, vector<int> & numSps, vector<int> & types)
{
    int sum=0;
    for(int i=0;i<prt->numOfsps.size();i++)
    {
        //if(prt->canos[i]==1 && prt->types[i]<2)
        //if(prt->canos[i]==1 && (prt->types[i]<2 || prt->types[i]>=3)) //Support splice variants
		if(prt->canos[i]==1) //Include read-throughs 
		{
            split_rna_t st=prt->sprnas[sum];
            st_vc.push_back(st);
            numSps.push_back(prt->numOfsps[i]);
			types.push_back(prt->types[i]);
        } else {
			split_rna_t st=prt->sprnas[sum];
			st_vc.push_back(st);
			numSps.push_back(prt->numOfsps[i]);
			types.push_back(-1);
		}
        sum+=prt->numOfsps[i];
    }
    return 0;
}


int Result::getAllJunctionsStep4(Gene& g, Reference & ref) {
    fjtvec2.clear();

    int index=0;
    for(int i=0;i<results.size();i++)
    {
        result_t * prt=this->getOneResult(i);
		
        if(prt->tier>3)
            break;
        if(prt->realPrint==0)
            continue;

		vector<int> types; //New

        vector<split_rna_t> st_vc;
        vector<int> numSps;
        bestSprnaVec2(prt,st_vc,numSps,types); //modified call
        //cout<<"Step4 size: "<<st_vc.size()<<endl;
		//for (int x=0;x<prt->types.size();x++)
		//	cout<<prt->types[x]<<endl;
		//cout<<""<<endl;
		//for (int x=0;x<types.size();x++)
		//	cout<<types[x]<<endl;
        for(int j=0;j<st_vc.size();j++)
        {
			//cout<<j<<endl;
			if(types[j]!=-1) { //==-1 if non-canonical
            	split_rna_t st=st_vc[j];
				int curr_type = types[j];
            	int p1;
            	if(st.bkLeft1==1)
                	p1=st.pos1;
            	else
                	p1=st.pos1+st.len1-1;


            	int p2;
            	if(st.bkLeft2==1)
             	   p2=st.pos2;
           	 	else
           	     p2=st.pos2+st.len2-1;

          		vector<junction_t> j1;
				vector<junction_t> j2;

				
          		g.getBestExon2(st.geneId1,p1, st.bkLeft1,j1);
          		g.getBestExon2(st.geneId2,p2, st.bkLeft2,j2); 
				


            	for (int x=0;x<j1.size(); x++)
            	{
                	for(int y=0;y<j2.size();y++)
                	{

                    	fusion_junction_2_t fjt;
						fjt.type = curr_type;
                    	fjt.fusion_id=index+1;
						
                    	fjt.tier=prt->tier;
                    	fjt.numSpReads=numSps[j];
						//fjt.numEnReads=prt->numOfEnRna;
                    	fjt.isCano=1;

                    	if(j1[x].is5p==1 && j2[y].is5p==0)
                    	{
                    	    getRefExonSeq(ref,j1[x].tid,j1[x].pos1+1,j1[x].pos2,j1[x].strand,fjt.seq1);
                    	    getRefExonSeq(ref,j2[y].tid,j2[y].pos1+1,j2[y].pos2,j2[y].strand,fjt.seq2);
                    	    fjt.p5=j1[x];
                    	    fjt.p3=j2[y];

                    	    if(fjt.p5.strand==1)
                    	        fjt.p5.pos1=fjt.p5.pos1+1;
                    	    if(fjt.p3.strand==0)
                    	        fjt.p3.pos1=fjt.p3.pos1+1;

                    	}
                    	else if(j1[x].is5p==0 && j2[y].is5p==1)
                    	{

                    	    getRefExonSeq(ref,j2[y].tid,j2[y].pos1+1,j2[y].pos2,j2[y].strand,fjt.seq1);
                    	    getRefExonSeq(ref,j1[x].tid,j1[x].pos1+1,j1[x].pos2,j1[x].strand,fjt.seq2);
                    	    fjt.p5=j2[y];
                    	    fjt.p3=j1[x];
                    	    if(fjt.p5.strand==1)
                    	        fjt.p5.pos1=fjt.p5.pos1+1;
                    	    if(fjt.p3.strand==0)
                    	        fjt.p3.pos1=fjt.p3.pos1+1;
                    	}
                    	fjtvec2.push_back(fjt);

                	}
            	}
			}
        }
        index++;
    }
    return 0;
}

int bestSprnaVec3(result_t * prt,vector<split_rna_t> & st_vc, vector<int> & numSps, vector<int> & types)
{
    int sum=0;
    for(int i=0;i<prt->numOfsps.size();i++)
    {
        //if(prt->canos[i]==0 && prt->types[i]<2)
        //if(prt->canos[i]==0 && (prt->types[i]<2 || prt->types[i]>=3))
		if(prt->canos[i]==0) //Include read-throughs
		{
            split_rna_t st=prt->sprnas[sum];
            st_vc.push_back(st);
            numSps.push_back(prt->numOfsps[i]);
			types.push_back(prt->types[i]);
        } else {
			split_rna_t st=prt->sprnas[sum];
			st_vc.push_back(st);
			numSps.push_back(prt->numOfsps[i]);
			types.push_back(-1);
		}
        sum+=prt->numOfsps[i];
    }
    return 0;
}

int Result::getAllJunctionsStep5(Gene& g, Reference & ref) {
    int index=0;
    for(int i=0;i<results.size();i++)
    {
        result_t * prt=this->getOneResult(i);
        if(prt->realPrint==0)
            continue;

		vector<int> types; //New

        vector<split_rna_t> st_vc;
        vector<int> numSps;
        bestSprnaVec3(prt,st_vc,numSps, types);
		//for (int x=0;x<prt->types.size();x++)
		//	cout<<prt->types[x]<<endl;
		//cout<<""<<endl;
		//for (int x=0;x<types.size();x++)
		//	cout<<types[x]<<endl;
        for(int j=0;j<st_vc.size();j++)
        {
			if(types[j]!=-1) {
            	split_rna_t st=st_vc[j];
				//int curr_type = types[j];
            	int p1;
            	if(st.bkLeft1==1)
                	p1=st.pos1;
            	else
                	p1=st.pos1+st.len1-1;


            	int p2;
            	if(st.bkLeft2==1)
                	p2=st.pos2;
            	else
                	p2=st.pos2+st.len2-1;

            	fusion_junction_2_t fjt;
            	fjt.fusion_id=index+1;
           
            	fjt.tier=prt->tier;
            	fjt.numSpReads=numSps[j];
				//fjt.numEnReads=prt->numOfEnRna;
            	fjt.isCano=0;
				fjt.type=types[j];
           
            	fjt.st=st;
            	if(g.isAt5p(st.geneId1,st.bkLeft1)==1)
            	{
                	fjt.geneId5p=st.geneId1;
                	fjt.geneId3p=st.geneId2;
                	fjt.pos5p=p1;
                	fjt.pos3p=p2;
            	}
            	else
            	{
                	fjt.geneId5p=st.geneId2;
                	fjt.geneId3p=st.geneId1;
                	fjt.pos5p=p2;
                	fjt.pos3p=p1;
            	}

            	fjtvec2.push_back(fjt);
			}
        }
        index++;
    }
    return 0;
}

bool my_sort_for_smc(fusion_junction_2_t i, fusion_junction_2_t j) // id, junc5p, junc3p,in-frame, short, short
{
	
    if(i.fusion_id<j.fusion_id)
    {
        return true;
    }
    else if(i.fusion_id==j.fusion_id)
    {
        int p5_pos1;
        if(i.isCano==1)
        {
            if(i.p5.strand==0)
                p5_pos1=i.p5.pos2;
            else
                p5_pos1=i.p5.pos1;
        }
        else
        {
            p5_pos1=i.pos5p;
        }

        int p5_pos2;
        if(j.isCano==1)
        {
            if(j.p5.strand==0)
                p5_pos2=j.p5.pos2;
            else
                p5_pos2=j.p5.pos1;
        }
        else
        {
            p5_pos2=j.pos5p;
        }

                if(p5_pos1<p5_pos2)
                {
                    return true;
                }
                else if(p5_pos1==p5_pos2)
                {

                    int p3_pos1;
                    if(i.isCano==1)
                    {
                        if(i.p3.strand==0)
                            p3_pos1=i.p3.pos1;
                        else
                            p3_pos1=i.p3.pos2;
                    }
                    else
                    {
                        p3_pos1=i.pos3p;
                    }

                    int p3_pos2;
                    if(j.isCano==1)
                    {
                    if(j.p3.strand==0)
                        p3_pos2=j.p3.pos1;
                    else
                        p3_pos2=j.p3.pos2;
                    }
                    else
                    {
                        p3_pos2=j.pos3p;
                    }    

                    if(p3_pos1<p3_pos2)
                    {
                        return true;
                    }
                    else if(p3_pos1==p3_pos2)
                    {
						if(i.numSpReads>j.numSpReads) { //New.
							return true;
						}
						else 
						{ //Original
                        	if(i.isInframe>j.isInframe)
                        	{
                            	return true;
                        	}
                        	else if(i.isInframe==j.isInframe)
                        	{
                            	if(i.p5.pos2-i.p5.pos1 < j.p5.pos2-j.p5.pos1)
                            	{
                                	return true;
                            	}
                            	else if(i.p5.pos2-i.p5.pos1 == j.p5.pos2-j.p5.pos1)
                            	{
                                	if(i.p3.pos2-i.p3.pos1 < j.p3.pos2-j.p3.pos1)
                                	{
                                    	return true;
                                	}
                                	else
                                	{
                                    	return false;
                                	}
                            	}
                            	else
                            	{
                                	return false;
                            	}
                        	}
                        	else
                        	{
                            	return false;
                        	}
						}
                    }
                    else
                    {
                        return false;
                    }

                }
        else
        {
            return false;
        }

    }
    else
    {
        return false;
    }

}
int Result::getAllJunctionsStep6(char * file,Gene& g, Reference & ref, char* dir) {
    Updator updator;
    sort(fjtvec2.begin(),fjtvec2.end(),my_sort_for_smc);
	std::string filename = createFilename(dir,file);
    ofstream outFile(filename);

    int last_pos5=-1;
    int last_pos3=-1;
	int last_support=-1;

	
    for (int i=0; i<fjtvec2.size(); i++)
    {
		//cout<<"i="<<i<<" isCano"<<fjtvec2[i].isCano<<endl;
        uint32_t p5_pos;
        if(fjtvec2[i].isCano==1)
        {
            if(fjtvec2[i].p5.strand==0)
                p5_pos=fjtvec2[i].p5.pos2;
            else
                p5_pos=fjtvec2[i].p5.pos1;
        }
        else
        {
            p5_pos=fjtvec2[i].pos5p;
        }

        uint32_t p3_pos;
        if(fjtvec2[i].isCano==1)
        {
        if(fjtvec2[i].p3.strand==0)
            p3_pos=fjtvec2[i].p3.pos1;
        else
            p3_pos=fjtvec2[i].p3.pos2;
        }
        else
        {
            p3_pos=fjtvec2[i].pos3p;
        }
		if(fjtvec2[i].isCano==0)
                updator.update(fjtvec2[i].geneId5p, fjtvec2[i].geneId3p, p5_pos, p3_pos, fjtvec2[i].st, ref, g);

        //New
		//Confirm pos is valid based on gene location
		
		if(fjtvec2[i].isCano==1) {
			if(!g.posIsValid(fjtvec2[i].p5.gId, p5_pos)) {
				int tmp = p5_pos;
				p5_pos = p3_pos;
				p3_pos = tmp;
			}
		} else {
			
			if(!g.posIsValid(fjtvec2[i].geneId5p, p5_pos)) {
				int tmp = p5_pos;
				p5_pos = p3_pos;
				p3_pos = tmp;
			}
		}
		
		if(last_pos5!=p5_pos || last_pos3!=p3_pos )
		//Edge case where two junctions are extremely close together and get mapped similarly
		//if(last_pos5!=p5_pos || last_pos3!=p3_pos || last_support != fjtvec2[i].numSpReads)
        {
            if(fjtvec2[i].isCano==1)
                outFile<<ref.getCharName(fjtvec2[i].p5.tid)<<"\t";
            else
                outFile<<ref.getCharName(g.getTid(fjtvec2[i].geneId5p))<<"\t";
            int strand1;
            if(fjtvec2[i].isCano==1)
                strand1=fjtvec2[i].p5.strand;
            else
                strand1=g.getStrand(fjtvec2[i].geneId5p);

            if(strand1==0)
            {
                outFile<<"-1"<<"\t";
                outFile<<p5_pos<<"\t";
            }
            else
            {
                outFile<<p5_pos-1<<"\t";
                outFile<<"-1"<<"\t";
            }
            if(fjtvec2[i].isCano==1)
                outFile<<ref.getCharName(fjtvec2[i].p3.tid)<<"\t";            
            else
                outFile<<ref.getCharName(g.getTid(fjtvec2[i].geneId3p))<<"\t";

            int strand2;
            if(fjtvec2[i].isCano==1)
                strand2=fjtvec2[i].p3.strand;
            else
                strand2=g.getStrand(fjtvec2[i].geneId3p); 
            if(strand2==0)
            {
                outFile<<p3_pos-1<<"\t";
                outFile<<"-1"<<"\t";
            }
            else
            {
                outFile<<"-1"<<"\t";
                outFile<<p3_pos<<"\t";
            }
            if(fjtvec2[i].isCano==1)
                outFile<<g.getName2(fjtvec2[i].p5.gId)<<">>"<<g.getName2(fjtvec2[i].p3.gId)<<"\t";
            else
                outFile<<g.getName2(fjtvec2[i].geneId5p)<<">>"<<g.getName2(fjtvec2[i].geneId3p)<<"\t";

            outFile<<fjtvec2[i].tier<<"\t";
            if(fjtvec2[i].isCano==1)
            {
                if(fjtvec2[i].p5.strand==0)
                    outFile<<"+\t";
                else
                    outFile<<"-\t";
            }
            else
            {
                int s=g.getStrand(fjtvec2[i].geneId5p);
                if(s==0)
                    outFile<<"+\t";
                else
                    outFile<<"-\t";
            }
            if(fjtvec2[i].isCano==1)
            {
                if(fjtvec2[i].p3.strand==0)
                    outFile<<"+\t";
                else
                    outFile<<"-\t";
            }
            else
            {
                int s=g.getStrand(fjtvec2[i].geneId3p);
                if(s==0)
                    outFile<<"+\t";
                else
                    outFile<<"-\t";
            }
            
			outFile<<fjtvec2[i].numSpReads<<"\t";
			outFile<<getType(fjtvec2[i].type)<<endl; //New column
			last_support = fjtvec2[i].numSpReads;

        }
        last_pos5=p5_pos;
        last_pos3=p3_pos;
    }



    return 0;


}


int Result::printFcirc(char* file, Gene & g, Reference & ref, char* dir) {
	Updator updator;


	vector<fcirc_t> fcircs;
	int last_pos5=-1;
	int last_pos3=-1;
	//Collect main junctions into fcircs vector
	for (int i=0;i<fjtvec2.size();i++) {
		//Done exactly as in Step6, to ensure consistency
		uint32_t p5_pos;
		if(fjtvec2[i].isCano==1) {
			if(fjtvec2[i].p5.strand==0)
				p5_pos=fjtvec2[i].p5.pos2;
			else
				p5_pos=fjtvec2[i].p5.pos1;
		} else {
			p5_pos=fjtvec2[i].pos5p;
		}
		
		uint32_t p3_pos;
		if(fjtvec2[i].isCano==1) {
			if(fjtvec2[i].p3.strand==0)
				p3_pos=fjtvec2[i].p3.pos1;
			else
				p3_pos=fjtvec2[i].p3.pos2;
		} else {
			p3_pos=fjtvec2[i].pos3p;
		}
		//Don't think this is necessary
		if(fjtvec2[i].isCano==0)
			updator.update(fjtvec2[i].geneId5p,fjtvec2[i].geneId3p,p5_pos,p3_pos,fjtvec2[i].st,ref,g);
		
		if(fjtvec2[i].isCano==1) {
			if(!g.posIsValid(fjtvec2[i].p5.gId,p5_pos)) {
				int tmp = p5_pos;
				p5_pos = p3_pos;
				p3_pos = tmp;
			}
		} else {
			if(!g.posIsValid(fjtvec2[i].geneId5p,p5_pos)) {
				int tmp = p5_pos;
				p5_pos = p3_pos;
				p3_pos = tmp;
			}
		}
		if(last_pos5!=p5_pos || last_pos3!=p3_pos) {
			//Only create potential fcircs for main junctions
			if(fjtvec2[i].type==0 || fjtvec2[i].type==1 || fjtvec2[i].type==2) {
				fcirc_t fusion;
				fusion.curr_id = fjtvec2[i].fusion_id;
				if(fjtvec2[i].isCano==1) {
					fusion.fusion5chrm = ref.getCharName(fjtvec2[i].p5.tid);
					fusion.fusion5strand = fjtvec2[i].p5.strand;
					fusion.fusion3chrm = ref.getCharName(fjtvec2[i].p3.tid);
					fusion.fusion3strand = fjtvec2[i].p3.strand;
					fusion.gene5 = g.getName2(fjtvec2[i].p5.gId);
					fusion.gene3 = g.getName2(fjtvec2[i].p3.gId);
				} else {
					fusion.fusion5chrm = ref.getCharName(g.getTid(fjtvec2[i].geneId5p));
					fusion.fusion5strand = g.getStrand(fjtvec2[i].geneId5p);
					fusion.fusion3chrm = ref.getCharName(g.getTid(fjtvec2[i].geneId3p));
					fusion.fusion3strand = g.getStrand(fjtvec2[i].geneId3p);
					fusion.gene5 = g.getName2(fjtvec2[i].geneId5p);
					fusion.gene3 = g.getName2(fjtvec2[i].geneId3p);
				}
				
				
				//fusion.mainNumSpReads = fjtvec2[i].numSpReads;
				fusion.fusion5breakPos = p5_pos;
				fusion.fusion3breakPos = p3_pos;
				fusion.name = fusion.gene5 + "--" + fusion.gene3;
				fusion.hasFcirc = 0;
				fcircs.push_back(fusion);
			}
		}
		last_pos5 = p5_pos;
		last_pos3 = p3_pos;
	} //End of creation of fcircs vector
	//Match backsplices to main junctions
	last_pos5 = -1;
	last_pos3 = -1;
	for (int i=0;i<fjtvec2.size();i++) {
		//For all backsplices
		if(fjtvec2[i].type==4) {
			//Repeat, as in Step6
			uint32_t p5_pos;
			if(fjtvec2[i].isCano==1) {
				if(fjtvec2[i].p5.strand==0)
					p5_pos=fjtvec2[i].p5.pos2;
				else
					p5_pos=fjtvec2[i].p5.pos1;
			} else {
				p5_pos=fjtvec2[i].pos5p;
			}
		
			uint32_t p3_pos;
			if(fjtvec2[i].isCano==1) {
				if(fjtvec2[i].p3.strand==0)
					p3_pos=fjtvec2[i].p3.pos1;
				else
					p3_pos=fjtvec2[i].p3.pos2;
			} else {
				p3_pos=fjtvec2[i].pos3p;
			}
			//Don't think this is necessary
			if(fjtvec2[i].isCano==0)
				updator.update(fjtvec2[i].geneId5p,fjtvec2[i].geneId3p,p5_pos,p3_pos,fjtvec2[i].st,ref,g);
		
			if(fjtvec2[i].isCano==1) {
				if(!g.posIsValid(fjtvec2[i].p5.gId,p5_pos)) {
					int tmp = p5_pos;
					p5_pos = p3_pos;
					p3_pos = tmp;
				}
			} else {
				if(!g.posIsValid(fjtvec2[i].geneId5p,p5_pos)) {
					int tmp = p5_pos;
					p5_pos = p3_pos;
					p3_pos = tmp;
				}
			}
			if(last_pos5!=p5_pos || last_pos3!=p3_pos) {
				string back5;
				string back3;
				if(fjtvec2[i].isCano==1) {
					back5 = g.getName2(fjtvec2[i].p5.gId);
					back3 = g.getName2(fjtvec2[i].p3.gId);
				} else {
					back5 = g.getName2(fjtvec2[i].geneId5p);
					back3 = g.getName2(fjtvec2[i].geneId3p);
				}
				for(int j=0;j<fcircs.size();j++){
					if(back5==fcircs[j].gene3 && back3==fcircs[j].gene5) {
						fcircs[j].fusion5backPos = p3_pos;
						fcircs[j].fusion3backPos = p5_pos;
						//fcircs[j].backspliceSpReads = fjtvec2[i].numSpReads;
						fcircs[j].hasFcirc = 1;
					} 
				}

			}
			last_pos5 = p5_pos;
			last_pos3 = p3_pos;
		}
	}
	//Print results
	std::string filename = createFilename(dir,file);
	ofstream outFile(filename);
	//Header
	outFile<<"#Id\tFusionName\tFusion5Prime.Backsplice\tFusion3Prime.Backsplice\tFusion5Prime.Breakpoint\tFusion3Prime.Breakpoint\n";
	for(int i=0;i<fcircs.size();i++) {
		if(fcircs[i].hasFcirc==1) {
			string strand5;
			string strand3;
			if(fcircs[i].fusion5strand==1)
				strand5="-";
			else
				strand5="+";
			if(fcircs[i].fusion3strand==1)
				strand3="-";
			else
				strand3="+";
			outFile<<fcircs[i].curr_id<<"\t";
			outFile<<fcircs[i].name<<"\t";
			outFile<<fcircs[i].fusion5chrm<<":"<<fcircs[i].fusion5backPos<<":"<<strand5<<"\t";
			outFile<<fcircs[i].fusion3chrm<<":"<<fcircs[i].fusion3backPos<<":"<<strand3<<"\t";
			outFile<<fcircs[i].fusion5chrm<<":"<<fcircs[i].fusion5breakPos<<":"<<strand5<<"\t";
			outFile<<fcircs[i].fusion3chrm<<":"<<fcircs[i].fusion3breakPos<<":"<<strand3<<"\n";
		}	
	}
	return 0;
}

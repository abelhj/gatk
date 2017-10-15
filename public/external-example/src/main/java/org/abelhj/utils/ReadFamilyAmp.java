package org.abelhj.utils;

import org.abelhj.utils.TypedTuple;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.ArrayList;
import java.util.List;
import java.util.LinkedList;
import java.util.Collections;
import java.util.Comparator;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;
import java.io.PrintStream;

public class ReadFamilyAmp {

    private String barcode=null;
    private String bestAmp=null;
    private Map<String , GATKSAMRecord[] > readfam=null;   //holds first and second read, indexed by name
    private int firstPos;                                  //start and chr of first occurence of this family
    private String firstChr;
    private String leadRead=null;                          //lead read name
    private boolean single;                                //family only contains one end (other end might be elsewhere, if improperly paired)
    private int single_order;                              //if single, order of lead read in pair

    public ReadFamilyAmp() {
	readfam=new LinkedHashMap<String, GATKSAMRecord[] >();
    }

    public ReadFamilyAmp(GATKSAMRecord rec) {
	this();
	bestAmp="None";
	barcode=rec.getStringAttribute("X0");
	firstPos=rec.getAlignmentStart();
	firstChr=rec.getReferenceName();
	add(rec);
	single=false;
	single_order=-1;
    }

    public String getLeadRead() {
	return leadRead;
    }

    public String getBestAmp() {
	return bestAmp;
    }

    public void add(GATKSAMRecord rec) {
	if(!readfam.containsKey(rec.getReadName())) {
	    readfam.put(rec.getReadName(), new GATKSAMRecord[2]);
	}
	if(rec.getFirstOfPairFlag()) {
	    readfam.get(rec.getReadName())[0]=rec;
	}
	else {
	    readfam.get(rec.getReadName())[1]=rec;
	}
    }
	
    public int getSize() {
	return readfam.size();
    }

    public String toString() {
	String str=barcode+"\t"+bestAmp+"\t[\n";
	for (String name : readfam.keySet()) {
	    for (int i=0; i<2; i++) {
		if(readfam.get(name)[i]!=null) {
		    str+=readfam.get(name)[i].getSAMString()+",\n";
		}
	    }
	}
	str+="]";
	return str;
    }

    //if both ends of any read not present, sets lead to null
    public boolean findLeadRead( ) {

	Map<Integer, Map<Integer, List<String> > > starts =new LinkedHashMap<Integer, Map<Integer, List<String> > >();
	boolean twoends=false;

	for (String rdname : readfam.keySet()) {
	    if(readfam.get(rdname)[0]!=null && readfam.get(rdname)[1]!=null) {
		twoends=true;
		int start0=readfam.get(rdname)[0].getAlignmentStart();
		int start1=readfam.get(rdname)[1].getAlignmentStart();
		if(!starts.containsKey(start0)) {
		    starts.put(start0, new LinkedHashMap<Integer, List<String> >());
		}
		if(!starts.get(start0).containsKey(start1)) {
		    starts.get(start0).put(start1, new ArrayList<String>());
		}
		starts.get(start0).get(start1).add(rdname);
	    }
	}

	if(!twoends) {
	    return twoends;
	}
	
	int maxct=-1;
	int freqstart0=-1;
	int freqstart1=-1;

	for(Integer start0: starts.keySet()) {
	    for(Integer start1 : starts.get(start0).keySet()) {
		int count=starts.get(start0).get(start1).size();
		if(count>maxct) {
		    maxct=count;
		    freqstart0=start0;
		    freqstart1=start1;
		}
	    }
	}

	Map<Integer, Map<Integer, List<String> > > lens=new LinkedHashMap<Integer, Map<Integer, List<String> > >();

	for(String rdname : starts.get(freqstart0).get(freqstart1)) {
	    int len0=readfam.get(rdname)[0].getReadLength();
	    int len1=readfam.get(rdname)[1].getReadLength();
	    
	    if(!lens.containsKey(len1)) {
		lens.put(len0, new LinkedHashMap<Integer, List<String> >());
	    }
	    if(!lens.get(len0).containsKey(len1)) {
		lens.get(len0).put(len1, new  ArrayList<String>());
	    }
	    lens.get(len0).get(len1).add(rdname);
	}

	maxct=-1;
	int freqlen0=-1;
	int freqlen1=-1;

	for(Integer len0: lens.keySet()) {
	    for(Integer len1 : lens.get(len0).keySet()) {
		int count=lens.get(len0).get(len1).size();
		if(count>maxct) {
		    maxct=count;
		    freqlen0=len0;
		    freqlen1=len1;
		}
	    }
	}

	int minscore=999;
	String bestRead="";

	for(String rdname : lens.get(freqlen0).get(freqlen1)) {
	    int curscore=readfam.get(rdname)[0].getIntegerAttribute("NM")+readfam.get(rdname)[1].getIntegerAttribute("NM");
	    if(curscore<minscore) {
		minscore=curscore;
		bestRead=rdname;
	    }
	}	
	if(minscore==999) {
	    return false;
	}
	leadRead=bestRead;
	return true;
    }

    public void findLeadReadSingle( Map<String, Integer> bcmaster) {

	boolean hasLead=false;
	for(String readname : readfam.keySet()) {
	    if(bcmaster.containsKey(readname)) {
		hasLead=true;
		leadRead=readname;
		bcmaster.remove(leadRead);
		for(int ii=0; ii<2; ii++) {
		    if(readfam.get(leadRead)[ii]!=null) {
			single_order=ii;
			break;
		    }
		}
		return;
	    }
	}

	int[] endcounts=new int[2];

	for(int ii=0; ii<2; ii++) {
	    for(String rdname : readfam.keySet()) {
		if(readfam.get(rdname)[ii]!=null) {
		    endcounts[ii]++;
		}
	    }
	}

	int pairEnd=0;
	if(endcounts[1]>endcounts[0]) {
	    pairEnd=1;
	}

	single_order=pairEnd;

        Map<Integer, List<String> >  starts=new LinkedHashMap<Integer, List<String> >  ();

	//System.err.println("single_order="+single_order);
	for (String rdname : readfam.keySet()) {
	    if(readfam.get(rdname)[single_order]!=null) {
		int start=readfam.get(rdname)[single_order].getAlignmentStart();
		//System.err.println(start+"\t"+rdname);
		if(!starts.containsKey(start)) {
		    starts.put(start, new ArrayList<String>());
		}
		starts.get(start).add(rdname);
	    }
	}

	int maxct=-1;
        int freqstart=-1;

        for(Integer start0: starts.keySet()) {
	    int count=starts.get(start0).size();
	    if(count>maxct) {
		maxct=count;
                freqstart=start0;
	    }
        }

	Map<Integer, List<String> > lens=new LinkedHashMap<Integer, List<String> > ();
	//System.err.println(toString());
	//System.err.println("freqstart="+freqstart);
	for(String rdname : starts.get(freqstart)) {
	    int len=readfam.get(rdname)[single_order].getReadLength();

	    if(!lens.containsKey(len)) {
		lens.put(len, new ArrayList<String>());
	    }
	    lens.get(len).add(rdname);
	}

	maxct=-1;
	int freqlen=-1;

	for(Integer len0: lens.keySet()) {
	    int count=lens.get(len0).size();
	    if(count>maxct) {
		maxct=count;
		freqlen=len0;
	    }
	}

	int minscore=999;
	String bestRead="";

	for(String rdname : lens.get(freqlen)) {
	    int curscore=readfam.get(rdname)[single_order].getIntegerAttribute("NM");
	    if(curscore<minscore) {
		minscore=curscore;
		bestRead=rdname;
	    }
	}		
	leadRead=bestRead;
	bcmaster.put(leadRead, 0);
    }


    public void findBestAmplicon(Map <String, Map<Integer, List< Amplicon> > > ampliconMap) {
	
	if(leadRead==null) {
	    System.err.println("Error: Must find lead read before best amplicon.\n");
	    System.exit(1);
	}

	int begin=-1;    //min and max aligned position for family 
	int end=-1;
	if(single) {
	    begin=readfam.get(leadRead)[single_order].getAlignmentStart();
	    end=readfam.get(leadRead)[single_order].getAlignmentEnd();
	} else {
	    //System.err.println(toString());
	    begin=readfam.get(leadRead)[0].getAlignmentStart();
	    end=readfam.get(leadRead)[1].getAlignmentEnd();
	    int begin1=readfam.get(leadRead)[1].getAlignmentStart();
	    int end1=readfam.get(leadRead)[0].getAlignmentEnd();
	    if(begin1<begin) {
		begin=begin1;
	    }
	    if(end1>end) {
		end=end1;
	    }
	}
	    
	String  curBest=null;
	int ampscore=1000;
	if(!ampliconMap.containsKey(firstChr)) {
	    return;   //amplicon is 'None'
	} else {
	    for(Integer endpos: ampliconMap.get(firstChr).keySet()) {
		if(Math.abs(endpos-end)<1000) {
		    for (Amplicon amp : ampliconMap.get(firstChr).get(endpos)) {
			int score=Math.abs(endpos-end)+Math.abs(amp.getStart()-begin);
			if(score<ampscore) {
			    ampscore=score;
			    curBest=amp.getName();;
			}
		    }
		}
	    }
	}
	bestAmp=curBest;
    }

    public boolean collapse(SAMFileWriter sw, Map <String, Map<Integer, List< Amplicon> >  > ampliconMap, 
			    Map<String, Integer> bcmaster ) {

	boolean paired=findLeadRead();
	if(paired) {
	    findBestAmplicon(ampliconMap);
	    System.out.println("collapsing\tlead="+leadRead+", bestAmplicon="+bestAmp);
	   
	    List< List<GATKSAMRecord> > recs=new ArrayList< List<GATKSAMRecord> >();

	    recs.add(new ArrayList<GATKSAMRecord>());
	    recs.add(new ArrayList<GATKSAMRecord>());
	    recs.get(0).add(readfam.get(leadRead)[0]);
	    recs.get(1).add(readfam.get(leadRead)[1]);
	    for(String nm : readfam.keySet()) {
		if(!nm.equals(leadRead)) {
		    if(readfam.get(nm)[0]!=null && readfam.get(nm)[0].getReadLength()==readfam.get(leadRead)[0].getReadLength()) {
			recs.get(0).add(readfam.get(nm)[0]);
		    }
		    if(readfam.get(nm)[1]!=null && readfam.get(nm)[1].getReadLength()==readfam.get(leadRead)[1].getReadLength()) {
			recs.get(1).add(readfam.get(nm)[1]);
		    }
		}
	    }
	    getConsensus(sw, recs.get(0));
	    getConsensus(sw, recs.get(1));
	    return true;
	} else {
	    single=true;
	    System.out.println("collapsing single end");
	    findLeadReadSingle(bcmaster);
	    findBestAmplicon(ampliconMap);
            System.out.println("collapsing\tlead="+leadRead+", bestAmplicon="+bestAmp);

	    List<GATKSAMRecord> recs=new ArrayList<GATKSAMRecord>();
	    recs.add(readfam.get(leadRead)[single_order]);
	    for(String nm : readfam.keySet()) {
		if(!nm.equals(leadRead)) {
		    if(readfam.get(nm)[single_order]!=null && readfam.get(nm)[single_order].getReadLength()==readfam.get(leadRead)[single_order].getReadLength()) {
			recs.add(readfam.get(nm)[single_order]);
		    }
		}
	    }
	    getConsensus(sw, recs);
			
	    return false;
	}
    }


    public void getConsensus( SAMFileWriter sw, List<GATKSAMRecord> reads) {
	
	StringBuilder bctag=new StringBuilder(""+readfam.size());
	boolean firstel=true;

	GATKSAMRecord lead=reads.get(0);
	int readlen=reads.get(0).getReadLength();
	byte[] newbases=new byte[readlen];
	byte[] newquals=new byte[readlen];

	for(int pos=0; pos<newbases.length; pos++) {
	    LinkedHashMap<Byte, ArrayList<TypedTuple<Integer, Byte> > > bases=new LinkedHashMap<Byte, ArrayList<TypedTuple<Integer, Byte> > >();
	    for(int rnum=0; rnum<reads.size(); rnum++) {
		GATKSAMRecord rr=reads.get(rnum);
		byte bb=rr.getReadBases()[pos];
		byte qq=rr.getBaseQualities()[pos];
		if(!bases.containsKey(bb)) {
		    bases.put(bb, new ArrayList<TypedTuple<Integer, Byte> >());
		}
		bases.get(bb).add(new TypedTuple<Integer, Byte>(rnum, qq));
	    }
	    int maxct=-1;
	    byte cbase=(byte)('N');
	    for(Byte bb : bases.keySet()) {
		if(bases.get(bb).size()>maxct) {
		    maxct=bases.get(bb).size();
		    cbase=bb;
		}
	    }
	    if(1.0*maxct/reads.size()<0.66) {
		newbases[pos]=(byte)('N');
		newquals[pos]=(byte)(5);
		cbase=(byte)('N');
	    } else {
		newbases[pos]=(byte)(cbase);
		int sum=0;
		int numc=0;
		for(TypedTuple<Integer, Byte> tt : bases.get(cbase)) {
		    sum+=tt.getRight().byteValue();
		    numc++;
		}
		newquals[pos]=QualityUtils.boundQual(sum/numc);
	    }
	    for(Byte bb: bases.keySet()) {
		if(bb!=cbase) {
		    for(TypedTuple<Integer, Byte> tt : bases.get(bb)) {
			if(firstel) {
			    bctag.append(":");
			    firstel=false;
			}
			bctag.append(pos+","+tt.getLeft()+","+(char)(bb.byteValue())+SAMUtils.phredToFastq(tt.getRight().byteValue())+";");
		    }
		}
	    }
	}
	lead.setBaseQualities(newquals);
	lead.setReadBases(newbases);
	lead.setAttribute("BC", bctag.toString());
	lead.setAttribute("X1", bestAmp);
	System.out.println(lead.getSAMString());
	sw.addAlignment(lead);
    }
}
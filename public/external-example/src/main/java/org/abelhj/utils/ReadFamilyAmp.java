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
    private Map<String , GATKSAMRecord[] > readfam=null;
    private int firstPos;
    private String firstChr;
    private String leadRead=null;

    public ReadFamilyAmp() {
	readfam=new LinkedHashMap<String, GATKSAMRecord[] >();
    }

    public ReadFamilyAmp(GATKSAMRecord rec) {
	this();
	bestAmp="None";
	barcode=rec.getStringAttribute("X0");
	firstPos=rec.getAlignmentStart();
	firstChr=rec.getReferenceName();
	readfam.put(rec.getReadName(), new GATKSAMRecord[2]);
	if(rec.getFirstOfPairFlag()) {
	    readfam.get(rec.getReadName())[0]=rec;
	}
	else {
	    readfam.get(rec.getReadName())[1]=rec;
	}
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
	String str=barcode+"\t[\n";
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

    public void findLeadRead( ) {

	int score=100;
	String best=null;

	for( String rdname : readfam.keySet()) {
	    if(readfam.get(rdname)[0]!=null && readfam.get(rdname)[1]!=null) {
		int nm1=readfam.get(rdname)[0].getIntegerAttribute("NM");
		int nm2=readfam.get(rdname)[1].getIntegerAttribute("NM");
		if(nm1+nm2<score) {
		    score=nm1+nm2;
		    best=rdname;
		}
	    }
	}
	leadRead=best;
    }


    public void findLeadReadSingle( ) {

	int score=100;
	String best=null;

	for( String rdname : readfam.keySet()) {
	    for(int i=0; i<2; i++) {
		if(readfam.get(rdname)[i]!=null ) {
		    int nm=readfam.get(rdname)[i].getIntegerAttribute("NM");
		    if(nm<score) {
			score=nm;
			best=rdname;
		    }
		}
	    }
	}
	leadRead=best;
    }


    public boolean findBestAmplicon(Map <String, Map<Integer, List< Amplicon> >  > ampliconMap) {
	
	if(leadRead==null) {
	    findLeadRead();
	    if(leadRead==null) {
		return false;
	    }
	}

	int begin=readfam.get(leadRead)[0].getAlignmentStart();
	int end=readfam.get(leadRead)[0].getAlignmentEnd();
	int begin1=readfam.get(leadRead)[1].getAlignmentStart();
	int end1=readfam.get(leadRead)[1].getAlignmentEnd();

	if(begin1<begin) {
	    begin=begin1;
	}
	if(end1>end) {
	    end=end1;
	}

	String  curBest=null;
	int ampscore=1000;
	if(!ampliconMap.containsKey(firstChr)) {
	    return true;
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
	return true;
    }

    public boolean collapse(SAMFileWriter sw, Map <String, Map<Integer, List< Amplicon> >  > ampliconMap, 
			    Map<String, Integer> bcmaster ) {

	findLeadRead();
	if(findBestAmplicon(ampliconMap)) {
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
	    System.out.println("collapsing single end");
	    boolean hasLead=false;
	    for(String readname : readfam.keySet()) {
		if(bcmaster.containsKey(readname)) {
		    hasLead=true;
		    leadRead=readname;
		    break;
		}
	    }
	    if(hasLead) {
		bcmaster.remove(leadRead);
	    } else {
		findLeadReadSingle();
		bcmaster.put(leadRead, 0);
	    }
	    int ii=0;
	    if(readfam.get(leadRead)[0]==null) {
		ii=1;
	    }
	    List<GATKSAMRecord> recs=new ArrayList<GATKSAMRecord>();
	    recs.add(readfam.get(leadRead)[ii]);
	    for(String nm : readfam.keySet()) {
		if(!nm.equals(leadRead)) {
		    if(readfam.get(nm)[ii]!=null && readfam.get(nm)[ii].getReadLength()==readfam.get(leadRead)[ii].getReadLength()) {
			recs.add(readfam.get(nm)[ii]);
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
	sw.addAlignment(lead);
    }
}
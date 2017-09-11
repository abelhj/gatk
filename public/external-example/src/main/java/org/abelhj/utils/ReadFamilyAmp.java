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
    private String baramp=null;
    private Map<String , List<GATKSAMRecord> > readfam=null;
    private int firstPos;
    private String firstChr;

    public ReadFamilyAmp() {
	readfam=new LinkedHashMap<String, List<GATKSAMRecord> >();
    }

    public ReadFamilyAmp(GATKSAMRecord rec) {
	this();
	barcode=read.getStringAttribute("X0");
	baramp=read.getStringAttribute("X1");
	firstPos=read.getAlignmentStart();
	firstChr=read.getReferenceName();
	readfam.put(rec.getReadName(), new ArrayList<GATKSAMRecord>(rec));
    }

    public void add(GATKSAMRecord rec) {
	if(!readfam.containsKey(rec.getName())) {
		readfam.put(rec.getName(), new ArrayList<GATKSAMRecord>(rec));
	} else {
	    readfam.get(rec.getName()).add(rec);
	}
    }
	

    public int getSize() {
	return readfam.size();
    }

    public String toString() {
	String str=barcode+"\t"+bamAmp+"\t[\n";
	for (String name : readfam.keySet()) {
	    for (GATKSAMRecord rec : readfam.get(name)) {
		str+=rec.getSAMString()+",\n";
	    }
	}
	str+="]";
	return str;
    }

    public String findLeadRead( Map<String, Map<Integer, Map<String, ReadFamilyLead> > >bcmaster) {
	String lead=null;
	if (bcmaster.containsKey(firstChr)) {
	    for(Integer ii : bcmaster.get(firstChr).keySet()) {
		if(abs(ii-firstPos)<500 && bcmaster.get(firstChr).get(ii).containsKey(barcode)) {
		    lead=bcmaster.get(firstChr).get(ii).get(barcode);
		    return lead;
		}
	    }
	}
	int score2=100;
	int score1=100;
	String best=null;
	boolean seentwo=false;
	for( String rdname : readfam.keySet()) {
	    if(readfam.get(rdname).size()==1) {
		int nm1=readfam.get(rdname).get(0).getStringAttribute("NM");
		if(nm1<score1) {
		    score1=nm1;
		    best=rdname;
		}
	    } else if (readfam.get(rdname).size()==2) {
		seentwo=true;
		int nm1=readfam.get(rdname).get(0).getStringAttribute("NM");
		int nm2=readfam.get(rdname).get(1).getStringAttribute("NM");
		if(nm1+nm2<score1+score2) {
		    score1=nm1;
		    score2=nm2;
		    best=rdname;
		}

	    }
	}
	if(!bcmaster.containsKey(firstChr)) {
	    bcmaster.put(firstChr, new LinkedHashMap<Integer, Map<String, String> >());
	}
	if(!bcmaster.get(firstChr).containsKey(firstPos)) {
	    bcmaster.get(firstChr).put(firstPos, new LinkedHashMap<String, String>());
	}
	Map<String, String> tempMap=new LinkedHashMap<String, String>();
	tempMap.put(barcode, best);
	bcmaster.get(firstChr).get(firstPos).put(tempMap);
	return best;
    }

    /*public static Comparator<GATKSAMRecord> SRComparator = new Comparator<GATKSAMRecord>() {         
	@Override         
	public int compare(GATKSAMRecord rec1, GATKSAMRecord rec2) {             
	    return rec1.getReadName().compareTo(rec2.getReadName());
	}     
	};*/       


    public void getConsensus(PrintStream bcout, Map<String, Map<String, String> > bcmaster, SAMFileWriter sw) {
	
	StringBuilder bctag=new StringBuilder(""+readfam.size());
	boolean firstel=true;

        String lead=findLeadRead(bcmaster);

	List< List<GATKSAMRecord> > recs=new ArrayList< List<GATKSAMRecord> >();

	//List<GATKSAMRecord> first=new ArrayList<GATKSAMRecord>();
	//List<GATKSAMRecord> second=new ArrayList<GATKSAMRecord>();

	int nends=readfam.get(lead).size();
	if(nends==1) {
	    recs.add(new ArrayList<GATKSAMRecord>());
	    recs.get(0).add(readfam.get(lead).get(0));
	    for(String nm : readfam.keySet()) {
		if(!nm.equals(lead)) {
		    recs.get(0).add(nm);
		}
	    }
	} else if (nends==2) {
	    recs.add(new ArrayList<GATKSAMRecord>());
	    recs.add(new ArrayList<GATKSAMRecord>());
	    if(readfam.get(lead).get(0).getFirstOfPairFlag()) {
		recs.get(0).add(readfam.get(lead).get(0));
		recs.get(1).add(readfam.get(lead).get(1));
		//first.add(readfam.get(lead).get(0));
                //second.add(readfam.get(lead).get(1));
	    } else {
		recs.get(0).add(readfam.get(lead).get(1));
		recs.get(1).add(readfam.get(lead).get(0));
		//first.add(readfam.get(lead).get(1));
                //second.add(readfam.get(lead).get(0));
	    }
	    for(String nm : readfam.keySet()) {
		if(!nm.equals(lead)) {
		    for(GATKSAMRecord rec : readfam.get(nm)) {
			if(rec.getFirstOfPairFlag()) {
			    recs.get(0).add(rec);
			} else {
			    recs.get(1).add(rec);
			}
		    }
		}
	    }
	} 

	for(int i=0; i<nends; i++) {
	    for(int j=0; j<recs.get(i).size(); j++) {
		rr=recs.get(i).get(j);
		System.out.println(baramp+"\t"+barcode+"\t"+rr.getSAMString());
	    }
	}

	/*for(int i=0; i<nends; i++) {


	    ArrayList<GATKSAMRecord> badlenReads=new ArrayList<GATKSAMRecord>();
	    readlen=recs.get(i).get(0).getReadLength();
	    for(int ii=0; ii<recs.get(i).size(); ii++) {
		rr=recs.get(i).get(ii);
		if(rr.getReadLength()!=readlen) {
		    badlenReads.add(rr);
		}
	    }
	}
	for(GATKSAMRecord rr : readfam) {
	    if(rr.getReadLength()!=readlen) {
		badlenReads.add(rr);
	    }
	}
	for(GATKSAMRecord rr : badlenReads) {
	    readfam.remove(rr);
	}
	      
	if(bcout !=null) {
	    for(GATKSAMRecord rr : readfam ) {
		bcout.println(barcode+"\t"+rr.getReadString()+"\t"+rr.getBaseQualityString());
	    }
	}

	byte[] newbases=new byte[lead.getReadLength()];
        byte[] newquals=new byte[lead.getReadLength()];

	for(int pos=0; pos<newbases.length; pos++) {
	    LinkedHashMap<Byte, ArrayList<TypedTuple<Integer, Byte> > > bases=new LinkedHashMap<Byte, ArrayList<TypedTuple<Integer, Byte> > >();
	    for(int rnum=0; rnum<readfam.size(); rnum++) {
		GATKSAMRecord rr=readfam.get(rnum);
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
	    if(1.0*maxct/readfam.size()<0.66) {
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
	sw.addAlignment(lead);
        if(bcout !=null) {
	    bcout.println(barcode+"\t"+lead.getReadString()+"\t"+lead.getBaseQualityString()+"\t"+lead.getStringAttribute("BC"));
	    } */  
    }
}
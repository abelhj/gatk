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

public class ReadFamily {

    private String barcode=null;
    private LinkedList<GATKSAMRecord> readfam=null;

    public ReadFamily() {
    }

    public ReadFamily(LinkedList<GATKSAMRecord> reads) {
	readfam=reads;
        barcode=readfam.get(0).getStringAttribute("X0");
    }

    public int getSize() {
	return readfam.size();
    }

    public String toString() {
	String str=barcode+"\t[\n";
	for (GATKSAMRecord read : readfam) {
	    str+=read.getSAMString()+",\n";
	}
	str+="]";
	return str;
    }

    public static Comparator<GATKSAMRecord> SRComparator = new Comparator<GATKSAMRecord>() {         
	@Override         
	public int compare(GATKSAMRecord rec1, GATKSAMRecord rec2) {             
	    return rec1.getReadName().compareTo(rec2.getReadName());
	}     
    };       

    public void getConsensus(PrintStream bcout, LinkedHashMap< String, LinkedHashMap<String, ReadFamilyLead> >  bcmaster, SAMFileWriter sw) {
	
	StringBuilder bctag=new StringBuilder(""+readfam.size());
	boolean firstel=true;

	Collections.sort(readfam, SRComparator);  //sort in order of readname
	boolean chg=false;
	GATKSAMRecord lead=readfam.get(0);

	//check to see if mate of one of the reads has already been chosen as lead. If so, make it lead; if not, first is alphabetical order is lead.
	if(bcmaster.containsKey(barcode)) {
	    for (GATKSAMRecord rr : readfam) {
		if(bcmaster.get(barcode).containsKey(rr.getReadName())) {
		    chg=true;
		    lead=rr;
		    bcmaster.get(barcode).get(rr.getReadName()).addRead(rr);
		}
	    }
	} else {
	    bcmaster.put(barcode, new LinkedHashMap<String, ReadFamilyLead>());
	    bcmaster.get(barcode).put(lead.getReadName(), new ReadFamilyLead(barcode, lead));
	}
	//if master contains one of the reads, and it's not at beginning of list, put it there
	if(chg && !readfam.get(0).equals(lead)) {
	    readfam.remove(lead);
 	    readfam.addFirst(lead);
	}
	//just remove reads whose lengths don't match lead
	ArrayList<GATKSAMRecord> badlenReads=new ArrayList<GATKSAMRecord>();
	int readlen=lead.getReadLength();
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
        }   
    }
}
package org.abelhj.utils;

import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

public class Amplicon {
   
    protected int start;
    protected int end;
    protected String name;
    protected String chr;

    public Amplicon(String chr, int start, int end, String nm) {
        this.chr=chr;
	this.start=start;
	this.end=end;
        this.name=nm;
    }

    public String toString() {
	String str=chr+":"+start+"-"+end+":"+name;
        return str;
    }

    public int getStart() {
	return start;
    }

    public int getEnd() {
	return end;
    }

    public String getChrom() {
	return chr;
    }

    public String getName() {
	return name;
    }

    public int getDist(GATKSAMRecord rec) {
	int dist=0;
	if(rec.getReferenceName().equals(chr)) {
	    dist=Math.abs(rec.getAlignmentStart()-start);
            int dist2=Math.abs(rec.getMateAlignmentStart()-start);
	    if(dist2<dist)
		dist=dist2;
	    return dist;
	} else {
	    return 9999;
	}
    }

    /*    public int getDist(GATKSAMRecord rec) {
	int dist1=0;
	int dist2=0;
	if(rec.getReferenceName().equals(chr)) {
	    if(rec.getFirstOfPairFlag()) {
		if(!rec.getReadNegativeStrandFlag()) {
		    dist1=Math.abs(rec.getAlignmentStart()-start);
		    //dist2=Math.abs(rec.getAlignmentStart()+Math.abs(rec.getInferredInsertSize())-end);
		} else {
		    //dist1=Math.abs(rec.getMateAlignmentStart()-start);
		  dist1=Math.abs(rec.getMateAlignmentStart()+Math.abs(rec.getInferredInsertSize())-end);
		}
	    } else {
		if(!rec.getReadNegativeStrandFlag()) {
		    dist1=Math.abs(rec.getAlignmentStart()-start);
		    //dist2=Math.abs(rec.getAlignmentStart()+Math.abs(rec.getInferredInsertSize())-end);
		} else {
		    //dist1=Math.abs(rec.getMateAlignmentStart()-start);
		  dist1=Math.abs(rec.getMateAlignmentStart()+Math.abs(rec.getInferredInsertSize())-end);
		}
	    }
	    return dist1+dist2;
	} else {
	    return 9999;
	}

	}*/


}
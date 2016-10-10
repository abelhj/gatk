package org.abelhj.utils;

import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;


public class ReadFamilyLead {

    private int read1start;
    private int read2start;
    private String barcode;

    private String read1name;
    private String read2name;

    public ReadFamilyLead(String bc) {
	barcode=bc;
	read1start=-1;
	read2start=-1;
	read1name=null;
	read2name=null;
    }

    public ReadFamilyLead(String bc, GATKSAMRecord rec) {
	this(bc);
	addRead(rec);
    }

    public String toString() {
	String str=barcode+"\t"+read1name+"\t"+read1start+"\t"+read2name+"\t"+read2start;
	return str;
    }

    public void addRead(GATKSAMRecord rec) {
	if(rec.getFirstOfPairFlag()) {
            read1start=rec.getAlignmentStart();
            read1name=rec.getReadName();
	} else {
            read2start=rec.getAlignmentStart();
            read2name=rec.getReadName();
	}
    }

    public boolean hasRead1(String name) {
	return (!(read1name==null));
    }

    public boolean hasRead2(String name) {
	return (!(read2name==null));
    }

}
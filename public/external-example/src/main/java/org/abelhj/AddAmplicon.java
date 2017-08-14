package org.abelhj;
//package org.broadinstitute.gatk.utils.examples;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.abelhj.utils.Amplicon;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFileWriterFactory;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;


@Reference(window=@Window(start=-1, stop=1))
public class AddAmplicon extends ReadWalker<Integer,Integer>  {


    @Argument(fullName = "outBam", shortName = "outBam", doc="BAM ouput", required=true)
    String outBam=null;
    @Argument(fullName = "ampliconBed", shortName = "amps", doc="amplicon bedfile", required=true)
    String ampliconBed = null;

    SAMFileWriter samwriter=null;
    SAMFileWriterFactory sf=null;
    LinkedHashMap <String, LinkedHashMap<Integer, ArrayList<Amplicon> >  > ampliconMap = null;

    public void initialize() {

        ampliconMap=new LinkedHashMap <String, LinkedHashMap<Integer, ArrayList< Amplicon> >  >();
        try {
    	BufferedReader br = new BufferedReader(new FileReader(ampliconBed));
    	String line;
    	while ((line = br.readLine()) != null) {
    	    String[] spl=line.split("\t");
    	    String chr=spl[0];
    	    int start=Integer.parseInt(spl[1]);
    	    Amplicon amp=new Amplicon(chr, start, Integer.parseInt(spl[2]), spl[3]);
    	    if(!ampliconMap.containsKey(chr)) {
    		ampliconMap.put(chr, new LinkedHashMap<Integer, ArrayList<Amplicon> >());
    	    }
    	    if(!ampliconMap.get(chr).containsKey(start)) {
    		ampliconMap.get(chr).put(start, new ArrayList<Amplicon> ());
    	    }
    	    ampliconMap.get(chr).get(start).add(amp);
    	}
        }
        catch (Exception e) {
    	e.printStackTrace();
        }
        sf=new SAMFileWriterFactory();
        samwriter=sf.makeBAMWriter(this.getToolkit().getSAMFileHeader(), true, new File(outBam), 5);
    }

    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
	
	if(ref!=null) {
	    if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case
	    
	    String chr=read.getReferenceName();
	    int start=read.getAlignmentStart();
	    int mindist=2000;
	    Amplicon bestAmp=null;
	    for(int i: ampliconMap.get(chr).keySet()) {
		if(Math.abs(i-start)<1000) {
		    for(Amplicon amp : ampliconMap.get(chr).get(i)) {
			int dist=amp.getDist(read);
			if(dist<mindist) {
			    mindist=dist;
			    bestAmp=amp;
			}
		    }
		}
	    }
            read.setAttribute("X1", bestAmp.getName());
	    samwriter.addAlignment(read);
	}
	return 1;
    }
   
    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return 0;
    }

    public void onTraversalDone(Integer result) {
	samwriter.close();
    }
}

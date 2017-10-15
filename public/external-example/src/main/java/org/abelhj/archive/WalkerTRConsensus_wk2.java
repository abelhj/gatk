package org.abelhj;

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
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.abelhj.utils.TypedTuple;
import htsjdk.samtools.SAMUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.abelhj.utils.ReadFamily;
import org.abelhj.utils.ReadFamilyLead;
import org.abelhj.utils.Amplicon;
import org.abelhj.utils.ReadFamilyAmp;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFileWriterFactory;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;


@Reference(window=@Window(start=-1, stop=1))
public class WalkerTRConsensus_wk2 extends ReadWalker<Integer,Integer>  {

    @Output
    PrintStream out;

    @Argument(fullName = "consensusBam", shortName = "cbam", doc="consensus BAM ouput", required=true)
    String consensusBam=null;
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    @Argument(fullName = "debug", shortName = "debug", doc= "1 to print read groups and consensus to bcfile", required=false)
    boolean debug=false;
    @Argument(fullName = "bcfile", shortName = "bc", doc= "barcode list", required=false)
    String bcfile=null;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
    int maxNM=99;
    @Argument(fullName = "ampliconBed", shortName = "ampliconBed", doc = "amplicon bed file", required=true)
	String ampliconBed=null;

    PrintStream bcout=null;
    SAMFileWriter samwriter=null;
    SAMFileWriterFactory sf=null;

    Map<String, Map<Integer, Map<String, ReadFamilyAmp> > >bcreads=null;
    Map <String, Map<Integer, List< Amplicon> >  > ampliconMap;


    String oldchr=null;
    String curchr=null;
    int curpos=-1;
    int oldpos=-1;
    String oldAmplicon=null;
    String curAmplicon=null;
    
    public void initialize() {

	ampliconMap=new LinkedHashMap <String, Map<Integer, List< Amplicon> >  >();
	try {
	    BufferedReader br = new BufferedReader(new FileReader(ampliconBed));
	    String line;
	    while ((line = br.readLine()) != null) {
		String[] spl=line.split("\t");
		String chr=spl[0];
		int start=Integer.parseInt(spl[1]);
		int end=Integer.parseInt(spl[2]);
		Amplicon amp=new Amplicon(chr, start, end, spl[3]);
		if(!ampliconMap.containsKey(chr)) {
		    ampliconMap.put(chr, new LinkedHashMap<Integer, List<Amplicon> >());
		}
		if(!ampliconMap.get(chr).containsKey(end)) {
		    ampliconMap.get(chr).put(end, new ArrayList<Amplicon> ());
		}
		ampliconMap.get(chr).get(end).add(amp);
	    }
	}
	catch (Exception e) {
	    e.printStackTrace();
	}
	if(debug) {
	    if(bcfile==null) {
		bcfile=consensusBam+".debug.txt";
	    }
	    try {
		bcout=new PrintStream(new File(bcfile));
	    } catch(Exception e) {
		System.err.println("Could not create debug output file.\n");
	    }
	}
	/*for (String chr : ampliconMap.keySet()) {
	    for (Integer endpos : ampliconMap.get(chr).keySet()) {
		for(Amplicon amp : ampliconMap.get(chr).get(endpos)) {
		    System.err.println(amp);
		}
	    }
	    }*/

	bcreads=new LinkedHashMap<String,  Map<Integer, Map<String, ReadFamilyAmp > > >();
	sf=new SAMFileWriterFactory();
	samwriter=sf.makeBAMWriter(this.getToolkit().getSAMFileHeader(), true, new File(consensusBam), 5);
    }

    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
	
	if(ref!=null) {
	    if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; 
	    curchr=read.getReferenceName();
	    curpos=read.getAlignmentStart();
	    if(oldchr==null) {
		oldchr=curchr;
		oldpos=curpos;
	    }
	    System.err.println(oldchr+"\t"+curchr+"\t"+curpos);
	    if(!curchr.equals(oldchr)) {
		if(ampliconMap.containsKey(oldchr)) {
		    System.err.println(oldchr+"\t"+curchr);
		    Map<Integer, List<String> > toremove=new LinkedHashMap<Integer, List<String> >();
		    for(Integer approxPos: bcreads.get(oldchr).keySet()) {
			for(String bc: bcreads.get(oldchr).get(approxPos).keySet()) {
			    System.out.println(curchr+"\t"+curpos+"\t"+approxPos+"\t"+bc);
			    ReadFamilyAmp rf=bcreads.get(oldchr).get(approxPos).get(bc);
			    System.out.println(rf);
			    if(rf.collapse(samwriter, ampliconMap)) {
				if(!toremove.containsKey(approxPos)) {
                                    toremove.put(approxPos, new ArrayList<String>());
                                }
                                toremove.get(approxPos).add(bc);
			    }
			}
		    }
		    for(Integer approxPos : toremove.keySet()) {
			for(String bc : toremove.get(approxPos))  {
			    bcreads.get(oldchr).get(approxPos).remove(bc);
			    if(bcreads.get(oldchr).get(approxPos).size()==0) {
				bcreads.get(oldchr).remove(approxPos);
			    }
			}
		    }
		    toremove=null;
		}
		oldchr=curchr;
	    } else if (bcreads.containsKey(curchr) && curpos!=oldpos) {
		Map<Integer, List<String> > toremove=new LinkedHashMap<Integer, List<String> >();
		for(Integer approxPos : bcreads.get(curchr).keySet()) {
		    if(curpos>approxPos*1000+2500) {
			for(String bc: bcreads.get(curchr).get(approxPos).keySet()){
			    System.out.println(curchr+"\t"+curpos+"\t"+approxPos+"\t"+bc);
			    ReadFamilyAmp rf=bcreads.get(curchr).get(approxPos).get(bc);
			    System.out.println(rf);
			    if(rf.collapse(samwriter, ampliconMap)){
				if(!toremove.containsKey(approxPos)) {
				    toremove.put(approxPos, new ArrayList<String>());
				}
				toremove.get(approxPos).add(bc);
			    }
			}   
			
		    }
		}
		for(Integer approxPos : toremove.keySet()) {
		    for(String bc : toremove.get(approxPos))  {
			bcreads.get(curchr).get(approxPos).remove(bc);
			if(bcreads.get(curchr).get(approxPos).size()==0) {
			    bcreads.get(curchr).remove(approxPos);
			}
		    }
		}
		toremove=null;
	    }
	    String bc=read.getStringAttribute("X0");
	    int approxPos=curpos/1000;
	    if(!bcreads.containsKey(curchr)) {
		bcreads.put(curchr, new LinkedHashMap<Integer, Map<String, ReadFamilyAmp> >());
	    }
      	    if(!bcreads.get(curchr).containsKey(approxPos) && bcreads.get(curchr).containsKey(approxPos-1)) {
		approxPos=approxPos-1;
	    } else if (!bcreads.get(curchr).containsKey(approxPos) && bcreads.get(curchr).containsKey(approxPos+1)) {
		approxPos=approxPos+1;
	    }
	    if(!bcreads.get(curchr).containsKey(approxPos)) {
		bcreads.get(curchr).put(approxPos, new LinkedHashMap<String, ReadFamilyAmp >());
	    }	
	    if(!bcreads.get(curchr).get(approxPos).containsKey(bc)) {
		ReadFamilyAmp rf = new ReadFamilyAmp(read);
		bcreads.get(curchr).get(approxPos).put(bc, rf);
	    } else {
		bcreads.get(curchr).get(approxPos).get(bc).add(read);
	    }
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
	/*for(String bc : bcreads.keySet()) {
	    for (String ii : bcreads.get(bc).keySet()) {
		ReadFamilyAmp rf=bcreads.get(bc).get(ii);
		//rf.getConsensus(bcout, bcmaster, samwriter);
	    }
	    }*/
	samwriter.close();
    }

}

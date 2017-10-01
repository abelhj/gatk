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
public class WalkerTRConsensus_wk1 extends ReadWalker<Integer,Integer>  {

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

    Map<String, Map<Integer, Map<String, ReadFamilyLead> > >bcmaster=null;
    Map<String, Map<String, ReadFamilyAmp >  >bcreads;
    Map <String, Map<Integer, List< Amplicon> >  > ampliconMap;


    String oldchr=null;
    String curchr=null;
    int oldpos=-1;
    int curpos=-1;
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
	bcreads=new LinkedHashMap<String,  Map<String, ReadFamilyAmp > >();
	bcmaster=new LinkedHashMap<String, Map<Integer, Map<String, ReadFamilyLead> > >();
	sf=new SAMFileWriterFactory();
	samwriter=sf.makeBAMWriter(this.getToolkit().getSAMFileHeader(), true, new File(consensusBam), 5);
    }

    public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
	
	if(ref!=null) {
	    if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; 
	    curchr=read.getReferenceName();
	    curpos=read.getAlignmentStart();
	    if(  !(curchr.equals(oldchr)) || ampliconMap.get(curchr).containsKey(curpos-10)) {         //find last possible current amplicon
		for(Integer endpos : ampliconMap.get(oldchr).keySet()) {
		    if (curpos>endpos+10 || !(curchr.equals(oldchr))) {
			for(Amplicon amp : ampliconMap.get(oldchr).get(endpos)) {
			    Map<String, ReadFamilyAmp> perAmpBC=bcreads.get(amp.getName());
			    perAmpBC.collapse(samwriter, ampliconMap, bcmaster);
			    bcreads.remove(amp.getName());
			}
		    }	
		}
		oldpos=curpos;
		oldchr=curchr;
		//oldAmplicon=curAmplicon;
	    }
	    String bc=read.getStringAttribute("X0");
	    String amp=read.getStringAttribute("X1");
	    if(!bcreads.containsKey(amp)) {
		bcreads.put(amp, new LinkedHashMap<String, ReadFamilyAmp>());
	    }
	    if(!bcreads.get(amp).containsKey(bc)) {
		ReadFamilyAmp rf = new ReadFamilyAmp(read);
		bcreads.get(amp).put(bc, rf);
	    } else {
		bcreads.get(amp).get(bc).add(read);
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
	for(String bc : bcreads.keySet()) {
	    for (String ii : bcreads.get(bc).keySet()) {
		ReadFamilyAmp rf=bcreads.get(bc).get(ii);
		//System.out.println(rf);
		rf.getConsensus(bcout, bcmaster, samwriter);
	    }
	}
	samwriter.close();
    }

    /*public Amplicon getNearestAmplicon(Map <String, List<GATKSAMRecord> > onebcmap, <String, Map<Integer, List< Amplicon> >  >ampliconMap) {




      }*/
}

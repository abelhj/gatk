package org.abelhj;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.genotyper.DiploidGenotype;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.File;


@Reference(window=@Window(start=-1, stop=1))
public class WalkerTR2017 extends LocusWalker<Integer,Integer>  {

    @Output
    PrintStream out;
    /**
     * Reads with mapping quality values lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    /**
     * Bases with quality scores lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth", required = false, minValue = 0, maxValue = Byte.MAX_VALUE)
    byte minBaseQuality = -1;
    @Argument(fullName = "minPercentRG", shortName = "mprg", doc = "Minimum percent mutant bases per read group", required = false, minValue = 0, maxValue=1)
    double minPercentRG=0.9;
    @Argument(fullName = "minCountPerBC", shortName = "minCtBC", doc= "Minimum number reads to accept bar code family", required=false)
    int minCountPerBC=1;
    @Argument(fullName = "debug", shortName = "debug", doc= "1 for verbose output", required=false)
    int debug=0;
    @Argument(fullName = "discardN", shortName = "discardN", doc= "alt is most common non-ref base excluding Ns", required=false)
    int discardN=0;
    @Argument(fullName = "bcfile", shortName = "bc", doc= "barcode list", required=false)
    String bcfile=null;
    @Argument(fullName = "writebc", shortName = "writebc", doc="0/1 to print list of barcodes", required=false)
    int writebc=0;
    @Argument(fullName = "minOffset", shortName = "minOffset", doc="min offset from either end of read", required=false)
    int minOffset=0;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
    int maxNM=99;

    boolean skipN=false;
    PrintStream bcout=null;
    
    public void initialize() {
	if(discardN==1) {
	    skipN=true;
	}
	if(writebc==1) {
	    if(bcfile!=null) {
		try {
		    bcout=new PrintStream(new File(bcfile));
		} catch(Exception e) {
		    System.err.println("barcode file not found\n");
		}
	     } else {
		System.err.println("Error: Must provide name for barcode list file.\n");
		System.exit(1);
	     }
	}
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

        ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality);
	String loc=pileup.getLocation().toString();
        String [] chrpos=loc.split(":");
        //LinkedHashMap<String, ArrayList<BaseFlag> > bcreads=new LinkedHashMap<String, ArrayList <BaseFlag> >();
        //BaseFlagMap bfmap=new BaseFlagMap();        
	LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<BaseFlag> > > bcreads_amp=new LinkedHashMap< Integer, LinkedHashMap<String, ArrayList <BaseFlag> > >(); 
	LinkedHashMap<Integer, BaseFlagMap> bfmap_amp=new LinkedHashMap<Integer, BaseFlagMap>();
        Character refbase=(char)ref.getBase();
        
        for(PileupElement p : pileup) {
	    GATKSAMRecord pread=p.getRead();
	    //check to see how softclips are handled with offset
	    if(pread.getIntegerAttribute("NM")<maxNM && p.getOffset()>=minOffset && p.getOffset()<=pread.getReadLength()-minOffset) {
		String rg=pread.getStringAttribute("X0");
		//int flag=pread.getFlags();
		int flag=0;
		if(pread.getFirstOfPairFlag()) {
		    flag+=64;
		}
		if(pread.getSecondOfPairFlag()) {
		    flag+=128;
		}
		if(pread.getReadNegativeStrandFlag()) {
		    flag+=16;
		} 
		char bb=(char)p.getBase();
		int ucstart=pread.getSoftStart();
		if(!bcreads_amp.containsKey(ucstart)) {
		    bcreads_amp.put(ucstart, new  LinkedHashMap<String, ArrayList<BaseFlag> >());
		    bfmap_amp.put(ucstart, new BaseFlagMap());
		}
		if(!bcreads_amp.get(ucstart).containsKey(rg)) {
		    bcreads_amp.get(ucstart).put(rg, new ArrayList<BaseFlag>());
		}
		bcreads_amp.get(ucstart).get(rg).add(new BaseFlag(bb, flag));
		bfmap_amp.get(ucstart).add(bb, flag);           
	    }
        }

	int minPerAmp=50;
	for(Integer ucstart: bcreads_amp.keySet()) {
	    char alt=bfmap_amp.get(ucstart).maxBase(refbase, false);
	    if(alt=='N') {
		alt='.';
	    }
		
	    out.print(chrpos[0]+"\t"+chrpos[1]+"\t"+ucstart+"\t"+refbase+"\t"+alt+"\t"+pileup.getBases().length+"\t"+bfmap_amp.get(ucstart).printSums(refbase)+"\t"+bfmap_amp.get(ucstart).printSums(alt)+"\t"+bfmap_amp.get(ucstart).printSums()+"\t");
	    out.print(bcreads_amp.get(ucstart).keySet().size()+"\t");
	    if(writebc==1) {
		bcout.print(chrpos[0]+"\t"+chrpos[1]+"\t"+refbase+"\t"+alt+"\t");
	    }
	    boolean comma=false;
	
	    BaseFlagMap rgbfmap=new BaseFlagMap();

	    for(String rg : bcreads_amp.get(ucstart).keySet()) {
		ArrayList<BaseFlag> reads=bcreads_amp.get(ucstart).get(rg);
		BaseFlagMap bfm=new BaseFlagMap();
		bfm.fill(reads);
		char mc_orig=bfm.maxBase();
		char mc=mc_orig;
		int flag=reads.get(0).flag;
		int sum=bfm.sum(mc);
		if(!BaseUtils.isRegularBase((byte)mc) || sum<1 || sum*1.0/bfm.sum()<minPercentRG || sum<minCountPerBC) {
		    mc='N';
		}
		rgbfmap.add(mc, flag);
		if(debug==1) {
		    bfm.print();
		}
		if(writebc==1) {
		    if(comma) {
			bcout.print(",");
		    }
		    bcout.print(rg+"_"+mc+"_"+mc_orig+"_"+sum+"_"+sum*1.0/bfm.sum());
		    comma=true;
		}
		if(mc_orig=='N') {
		    bfm.print();
		}
	    }
	    if(writebc==1) {
		bcout.println();
	    }
	    char bcalt=rgbfmap.maxBase(refbase, skipN);
	    if(bcalt=='N') {
		bcalt='.';
	    }
	    out.print(bcalt+"\t"+rgbfmap.printSums(refbase)+"\t"+rgbfmap.printSums(bcalt)+"\t"+rgbfmap.printSums()+"\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");
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
    }
}

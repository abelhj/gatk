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
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.BaseUtils;

import java.io.BufferedReader;
import java.io.FileReader;
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

import org.abelhj.utils.BaseFlag1;
import org.abelhj.utils.BaseFlagMap1;

@Reference(window=@Window(start=-1, stop=1))
public class AminoAcidTyper extends LocusWalker<Integer,Integer>  {

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
    @Argument(fullName = "debug", shortName = "debug", doc= "1 for verbose output", required=false)
    int debug=0;
    @Argument(fullName = "maxNM", shortName = "maxNM", doc="filter reads with edit distance greater than maxNM", required=false)
	int maxNM=99;
    @Argument(fullName = "invcf", shortName = "invcf", doc="per RG VAFs will be calculated only at variant positions", required=true)
	String invcf=null;

    LinkedHashMap<GenomeLoc, String> snps=null;
    GenomeLocParser gpl=null;
    ArrayList<Character> nucs=null;
    ArrayList<String> rgorder=null;

    public void initialize() {
        BufferedReader br=null;
        snps=new LinkedHashMap<GenomeLoc, String>();
        gpl=new GenomeLocParser(this.getMasterSequenceDictionary());
        try {
            br=new BufferedReader(new FileReader(invcf));
            String curLine=null;
            while((curLine=br.readLine())!=null) {
                String[] spl=curLine.split("\\t");
                GenomeLoc gl1=gpl.parseGenomeLoc(spl[0]+":"+spl[1]);
                snps.put(gl1, curLine);
            }
        }
        catch(Exception e) {
            e.printStackTrace();
        }
  	nucs=new ArrayList<Character>();
        nucs.add('A');
        nucs.add('C');
        nucs.add('G');
        nucs.add('T');
        nucs.add('N');
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case
        ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality);
	GenomeLoc loc=pileup.getLocation();
        if(!snps.containsKey(loc)) {
	    return null;
	} else {
   
	    LinkedHashMap<String, ArrayList<BaseFlag> > bcreads=new LinkedHashMap<String, ArrayList <BaseFlag> >();
	    BaseFlagMap bfmap=new BaseFlagMap();
	    //ArrayList<String> bcorder=null;        
	    Character refbase=(char)ref.getBase();
        
	    for(PileupElement p : pileup) {
		if(p.getRead().getIntegerAttribute("NM")<maxNM) {
		    String rg=p.getRead().getReadGroup().getId();
		    int flag=p.getRead().getFlags();
		    char bb=(char)p.getBase();
		    if(!bcreads.containsKey(rg)) {
			bcreads.put(rg, new ArrayList<BaseFlag>());
		    }
		    bcreads.get(rg).add(new BaseFlag(bb, flag));
		    bfmap.add(bb, flag);           
		}
	    }
        
	    char alt=bfmap.maxBase(refbase, false);
	    if(alt=='N') {
		alt='.';
	    }
		
	    out.print(loc+"\t"+refbase+"\t"+alt+"\t"+pileup.getBases().length+"\t"+bfmap.printSums(refbase)+"\t"+bfmap.printSums(alt)+"\t"+bfmap.printSums()+"\n");
	    if(rgorder==null) {
		rgorder=new ArrayList<String>(bcreads.keySet());
	    }

	    for(String rg : rgorder) {
		ArrayList<BaseFlag> reads=bcreads.get(rg);
		System.out.println("\t"+rg+"\t");
		for(BaseFlag bbff : bcreads.get(rg)) {
		    System.out.println(bbff);
		}
		System.out.println("\n\n");
		BaseFlagMap bfm=new BaseFlagMap();
		bfm.fill(reads);
		bfm.print();
		bfm.printSums();
		if(debug==1) {
		    bfm.print();
		}
	    }
	    /*boolean skipN=false;
	    char bcalt=rgbfmap.maxBase(refbase, skipN);
	    if(bcalt=='N') {
		bcalt='.';
	    }
	    out.print(bcalt+"\t"+rgbfmap.printSums(refbase)+"\t"+rgbfmap.printSums(bcalt)+"\t"+rgbfmap.printSums()+"\t"+(char)ref.getBases()[0]+"\t"+(char)ref.getBases()[2]+"\n");
	    */
	    return 1;
	}
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

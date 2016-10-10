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
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.BaseUtils;


import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

@Reference(window=@Window(start=-1, stop=1))
public class Haplotect extends LocusWalker<LocusNameAllele,Integer>  {

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
    @Argument(fullName = "pairsFile", shortName = "htp", doc = "file of locus pairs", required = true)
	String pairsFile;
    
    @Argument(fullName="debug", shortName = "dbg", doc="verbose for debugging", required=false)
    int debug=0;
    
    ArrayList<SnpPair> pairs=null;
    LinkedHashSet<GenomeLoc> snps=null;
    
    LinkedHashMap<SnpPair, LinkedHashMap<String, LinkedHashMap<Integer, PileupElement> > > hapmap=null;
    GenomeLocParser gpl=null;
    
    
    public void initialize() {
        BufferedReader br=null;
        pairs=new ArrayList<SnpPair>();
        snps=new LinkedHashSet<GenomeLoc>();
        gpl=new GenomeLocParser(this.getMasterSequenceDictionary());
        try {
            br=new BufferedReader(new FileReader(pairsFile));
            String curLine=null;
            while((curLine=br.readLine())!=null) {
                
                String[] spl=curLine.split("\\t");
                GenomeLoc gl1=gpl.parseGenomeLoc(spl[0]+":"+spl[1]);
                GenomeLoc gl2=gpl.parseGenomeLoc(spl[0]+":"+spl[2]);
                snps.add(gl1);
                snps.add(gl2);
                pairs.add(new SnpPair(gl1, gl2, curLine));
            }
        }
        catch(IOException e) {
            e.printStackTrace();
        }
        System.err.println("snp pairs");
        for(SnpPair sp:pairs) {
            System.err.println(sp);
        }
    }

    public LocusNameAllele map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

        ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality);
	    GenomeLoc loc=pileup.getLocation();
     
        if(snps.contains(loc)) {
            //System.err.println(loc);

	
            LinkedHashMap<String, PileupElement> readallele=new LinkedHashMap<String, PileupElement>();
            for(PileupElement p : pileup) {
                if(BaseUtils.isRegularBase(p.getBase())) 
                    readallele.put(p.getRead().getReadName(), p);
            }
            return new LocusNameAllele(loc, readallele);
        } else {
            return null;
        }

    }


    /**
     * 
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        
        hapmap=new LinkedHashMap<SnpPair, LinkedHashMap<String, LinkedHashMap<Integer, PileupElement> > >();
        for(SnpPair pr : pairs) 
            hapmap.put(pr, new LinkedHashMap<String, LinkedHashMap<Integer, PileupElement> >());
        return 0;
    }

    /**
     * 
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(LocusNameAllele value, Integer sum) {
        
        if(value!=null) {
       
            GenomeLoc loc=value.getLocus();
            for(SnpPair pr: pairs) {
                if(pr.matches(loc)) {
                    int order=pr.whichSnp(loc);
                    LinkedHashMap<String, PileupElement> nts=value.getAlleles();
                    for(String readname : nts.keySet()) {
                        if(!hapmap.get(pr).containsKey(readname)) 
                            hapmap.get(pr).put(readname, new LinkedHashMap<Integer, PileupElement>());
                        hapmap.get(pr).get(readname).put(order, nts.get(readname));
                    }    
                }
            }
        }
        return 0;
    }

 
    /**
     * when we finish traversing, close the result list
     * @param result the final reduce result
     */
    public void onTraversalDone(Integer result) {
        
        double aveContFrac=0;
        int totalCounts=0;
        int contamCounts=0;
        int informativeSites=0;
        LinkedHashMap<SnpPair, HapCounter> paircts=new LinkedHashMap<SnpPair, HapCounter>();
        for(SnpPair pr: hapmap.keySet()) {
            HapCounter hapcts=new HapCounter(pr);
            for(String name:hapmap.get(pr).keySet()) {
                if(hapmap.get(pr).get(name).keySet().size()==2) {               //only use reads covering both SNVs in pair
                    hapcts.add(hapmap.get(pr).get(name).get(1), hapmap.get(pr).get(name).get(2));
                }
            }
            paircts.put(pr, hapcts);
            System.out.print(pr.pairInfo()+"\t"+pr.distance()+"\t"+hapcts.totalCount()+"\t"+hapcts.countString());
            if(hapcts.numUniqueObsHaps()>2) {
                informativeSites++;
                totalCounts+=hapcts.totalCount();
                contamCounts+=hapcts.thirdHapCount();
                System.out.println("**");
               // hapcts.print();
            }
            else 
                System.out.println("  ");
            if(debug==1) {
                hapcts.print();
            }
                
        }
        aveContFrac=2.0*contamCounts/totalCounts;
        System.err.println("\n\nAve Contamination Fraction="+aveContFrac+"\tNum Informative Sites="+informativeSites);
      
        ArrayList<Double> alpha=new ArrayList<Double>();
        for(int i=0; i<5000; i++) {
            alpha.add(i/100000.0);
        }
        for(int i=6; i<50; i++) {
            alpha.add(i/100.0);
        }
        for(Double aa: alpha) {
            double loglik=0;
            for(SnpPair pr: paircts.keySet()) {
                HapCounter hapcts=paircts.get(pr);
                if(hapcts.totalCount()>0) {
                    loglik+=hapcts.getLogLik(aa);
                }
            }
            System.out.println(aa+"\t"+loglik);
        }
    }
}

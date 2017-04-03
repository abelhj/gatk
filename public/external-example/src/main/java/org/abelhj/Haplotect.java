package org.abelhj;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
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
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;
import java.lang.Runtime;

import org.abelhj.haplotect_utils.*;

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-1, stop=1))
public class Haplotect extends LocusWalker<LocusNameAllele,Integer>  {

    String version="0.2";
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
    @Argument(fullName = "debug", shortName = "dbg", doc = "verbose for debugging", required = false)
	int debug=0;
    boolean debugbool=false;
    
    @Argument(fullName = "gstol", shortName= "gstol", doc = "error tolerance for golden section search", required=false, minValue=0.00000001, maxValue=0.1)
	double gstol=0.005;
    @Argument(fullName = "outPrefix", shortName= "outPrefix", doc = "prefix for output files", required=true)
	String outPrefix=null;
	boolean uniform=false;
    @Argument(fullName = "unifPrior", shortName= "unif", doc = "use uniform prior for popn hap frequencies", required=false)
	int unifPrior=0;
    @Argument(fullName = "minreads", shortName = "mr", doc = "minium number of reads at a locus to include", required = false)
        int minreads=20;
    @Argument(fullName="minMultiPct", shortName = "minPct", doc="min pct multihaplotype sites", required=false)
	int minMultiPct=5;
    
    ArrayList<SnpPair> pairs=null;
    LinkedHashSet<GenomeLoc> snps=null;    
    LinkedHashMap<SnpPair, LinkedHashMap<String, LinkedHashMap<Integer, BaseandQual> > > hapmap=null;
    LinkedHashMap<SnpPair, HapCounter> paircts=null;
    HashSet<SnpPair> multiPairs=null;

    GenomeLocParser gpl=null;
    PrintStream out=null;
    PrintStream log=null;
    String currentContig=null;
    int totalCounts=0;
    int contamCounts=0;
    int informativeSites=0;
    int totalSites=0;
    double aveContFrac=0;
    double frac=0;
    double meancov=0;

    
    
    public void initialize() {
	if(unifPrior==1) {
	    uniform=true;
	}
	if(debug==1) {
	    debugbool=true;
	}
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
		//pairs.add(new SnpPair(gl1, gl2, curLine));
                pairs.add(new SnpPair(gl1, gl2, curLine, uniform));  //fix this !!
            }
	    br.close();
	    out=new PrintStream (new File(outPrefix+".txt"));
	    log=new PrintStream (new File(outPrefix+".multihaploci.txt"));
	    log.println("#Haplotect_v"+version); 
	    log.println("#"+pairs.size()+" SNP pairs read");
	    log.println("#chr\tSNP1\tSNP2\tall11\tall12\tall21\tall22\tpopn_counts\tdistance\ttotal_count\tsample_counts\tcontam_frac");
        }
        catch(IOException e) {
            e.printStackTrace();
        }
	//System.err.println("initialized\t"+Runtime.getRuntime().totalMemory()/(1024*1024)+"\t"+Runtime.getRuntime().freeMemory()/(1024*1024)+"\t"+Runtime.getRuntime().maxMemory()/(1024*1024));
	currentContig="";
    }

    public LocusNameAllele map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case
	GenomeLoc loc=context.getLocation();
	if(currentContig==null) {
	    currentContig=context.getContig();
        } else if(!context.getContig().equals(currentContig)) {
	    //System.err.println(loc+"\t"+Runtime.getRuntime().totalMemory()/(1024*1024)+"\t"+Runtime.getRuntime().freeMemory()/(1024*1024)+"\t"+Runtime.getRuntime().maxMemory()/(1024*1024));
            finishChrom();
	    currentContig=context.getContig();
	}
        if(snps.contains(loc)) {

	    ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(minBaseQuality, minMappingQuality);
            LinkedHashMap<String, BaseandQual> readallele=new LinkedHashMap<String, BaseandQual>();
            for(PileupElement p : pileup) {
                if(BaseUtils.isRegularBase(p.getBase())) 
                    readallele.put(p.getRead().getReadName(), new BaseandQual(p));
            }
            return new LocusNameAllele(loc, readallele);
        } else {
            return null;
        }
    }

    public Integer reduceInit() {
        
        hapmap=new LinkedHashMap<SnpPair, LinkedHashMap<String, LinkedHashMap<Integer, BaseandQual> > >();
	paircts=new LinkedHashMap<SnpPair, HapCounter>();
        multiPairs=new HashSet<SnpPair>();
        for(SnpPair pr : pairs) {
            hapmap.put(pr, new LinkedHashMap<String, LinkedHashMap<Integer, BaseandQual> >());
	}
	totalCounts=0;
        contamCounts=0;
        informativeSites=0;
        totalSites=0;
        aveContFrac=0;
        frac=0;
        meancov=0;
        return 0;
    }

    public Integer reduce(LocusNameAllele value, Integer sum) {
        
        if(value!=null) {
       
            GenomeLoc loc=value.getLocus();
            for(SnpPair pr: pairs) {
                if(pr.matches(loc)) {
                    int order=pr.whichSnp(loc);
                    LinkedHashMap<String, BaseandQual> nts=value.getAlleles();
                    for(String readname : nts.keySet()) {
                        if(!hapmap.get(pr).containsKey(readname)) 
                            hapmap.get(pr).put(readname, new LinkedHashMap<Integer, BaseandQual>());
                        hapmap.get(pr).get(readname).put(order, nts.get(readname));
                    }    
                }
            }
        }
        return 0;
    }

  public void finishChrom() {
      if(currentContig==null) {
	  return;
      }
    for(SnpPair pr: hapmap.keySet()) {
      HapCounter hapcts=new HapCounter(pr);
      for(String name:hapmap.get(pr).keySet()) {
        if(hapmap.get(pr).get(name).keySet().size()==2) {               //only use reads covering both SNVs in pair
          hapcts.add(hapmap.get(pr).get(name).get(1), hapmap.get(pr).get(name).get(2));
        }
      }
	    
      // requires minreads reads covering the locus to consider
      if (hapcts.totalCount()>minreads){		
        totalSites++;		
	paircts.put(pr, hapcts);
	meancov+=hapcts.totalCount();
	if(hapcts.numUniqueObsHaps()>2) {
	  double cFrac=0;
	  informativeSites++;	
	  totalCounts+=hapcts.totalCount();
	  multiPairs.add(pr);
	  if(hapcts.thirdHapCount()>0 && hapcts.fourthHapCount()==0){
	    contamCounts+=hapcts.thirdHapCount();
	    cFrac=2.0*hapcts.thirdHapCount()/hapcts.totalCount();
	    frac+=cFrac;
	  } else if(hapcts.fourthHapCount()>0) {
	    contamCounts+=hapcts.fourthHapCount();
	    cFrac=1.0*(hapcts.thirdHapCount()+hapcts.fourthHapCount())/hapcts.totalCount();
	    frac+=cFrac;
	  }		    
	  log.println(pr.pairInfo()+"\t"+pr.distance()+"\t"+hapcts.totalCount()+"\t"+hapcts.countString()+"\t"+String.format("%.4f",cFrac));
	} 
      }
      if(debug==1) {
        hapcts.print(log);
      }
    }
  }

  public void onTraversalDone(Integer result) {

	meancov/=totalSites;
	aveContFrac = 2.0*contamCounts/totalCounts;
	frac = frac / (1.0 * informativeSites);
	String[] spl=outPrefix.split("/");
	String id=spl[spl.length-1];
	PairLogLik mleCalculator=new PairLogLik(paircts, gstol, debugbool, log);
	mleCalculator.calcMleCI();	
	double mle=mleCalculator.getMle();
	double[] CI=mleCalculator.getCI();
	out.println("#sample\tmle\tCI_95\tmean_Frac\tmle_multi\tCI_95_multi\tNum_informative_sites\ttotalSites\tmean_Cov");
	if(1.0*informativeSites/totalSites>0.01*minMultiPct) {
	    mleCalculator=new PairLogLik(paircts, multiPairs, gstol, debugbool, log);
	    paircts=null;
	    mleCalculator.calcMleCI();
	    double mle_contam=mleCalculator.getMle();
	    double[] CI_contam=mleCalculator.getCI();
	    out.println(id+"\t"+String.format("%.3f",mle)+"\t"+String.format("%.3f", CI[0])+"-"+String.format("%.3f", CI[1])+"\t"+String.format("%.3f", aveContFrac)+"\t"+String.format("%.3f",mle_contam)+"\t"+String.format("%.3f", CI_contam[0])+"-"+String.format("%.3f", CI_contam[1])+"\t"+informativeSites+"\t"+totalSites+"\t"+String.format("%.3f",meancov));
	} else {
	    out.println(id+"\t"+String.format("%.3f",mle)+"\t"+String.format("%.3f", CI[0])+"-"+String.format("%.3f", CI[1])+"\t"+String.format("%.3f", aveContFrac)+"\tNA\tNA\t"+informativeSites+"\t"+totalSites+"\t"+String.format("%.3f",meancov));
	}

    }

}

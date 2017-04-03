package org.abelhj.haplotect_utils;

import java.util.LinkedHashMap;
import java.io.PrintStream;
import java.util.HashSet;

public class PairLogLik {

    private LinkedHashMap<SnpPair, HapCounter> paircts=null;
    private boolean debug;
    private double gstol=1;
    private double mle=-1;
    private double fmle=-1;
    private double cthresh=-1;
    private double[] CI={-1,-1};
    private double maxits=50;								//max number iterations for convergence of golden section search or bisection
    private double gold=(-1+Math.sqrt(5))/2;
    private LinkedHashMap<Double, Double> ff=null;			//stores calculated log-lik values
    boolean precise=true;
    private double ff_0=-99;
    private double ff_onehalf=-99;
    private PrintStream log=null;

    public PairLogLik(){}
	
    public PairLogLik( LinkedHashMap<SnpPair, HapCounter> paircts, double gstol, boolean debug, PrintStream log) {
	this.paircts=paircts;
	this.gstol=gstol;
	this.debug=debug;
	this.log=log;
	ff=new LinkedHashMap<Double, Double>();
    }

    public PairLogLik(LinkedHashMap<SnpPair, HapCounter> paircts1, HashSet<SnpPair> keepPairs, double gstol, boolean debug, PrintStream log) {
	paircts=new LinkedHashMap<SnpPair, HapCounter>();
	for(SnpPair sp: paircts1.keySet()) {
	    if(keepPairs.contains(sp)) {
		paircts.put(sp, paircts1.get(sp));
		//System.out.println(sp+"\t"+paircts1.get(sp).countString());
	    }
	}
        this.gstol=gstol;
        this.debug=debug;
        this.log=log;
        ff=new LinkedHashMap<Double, Double>();
    }
	
	public void calcMleCI() {
		calcMle();
		calcCI();
	}
	


	public double getMle() {
		return mle;
	}
	
	public double[] getCI() {
		return CI;
	}

	
	private void calcMle() {
	
		mle=0;
		fmle=0;
	
		//starting values for golden section search
		double aa=0;
		double bb=0.5;
		double xx=aa+gold*(bb-aa);
		double uu=(bb-aa)-xx;
		
		double faa=calcLL(aa);
		double fbb=calcLL(bb);
		ff_0=faa;
		ff_onehalf=fbb;
		double fxx=calcLL(xx);
		double fuu=calcLL(uu);
		
		ff.put(aa, faa);
		ff.put(bb, fbb);
		ff.put(xx, fxx);
		ff.put(uu, fuu);
		
		if(debug) {
			log.println("**0\t"+aa+"\t"+uu+"\t"+xx+"\t"+bb+"\t"+fuu+"\t"+fxx);
		}
		
		for(int i=1; i<=maxits; i++) {
			if(fuu>fxx) {
				mle=uu;
				fmle=fuu;
				bb=xx;
				xx=uu;
				fbb=fxx;
				fxx=fuu;
				uu=aa+gold*(xx-aa);
				fuu=calcLL(uu);
				ff.put(uu, fuu);
			} else {
				mle=xx;
				fmle=fxx;
				aa=uu;
				uu=xx;
				faa=fuu;
				fuu=fxx;
				xx=uu+gold*(bb-uu);
				fxx=calcLL(xx);
				ff.put(xx, fxx);
			}
			if(debug) {
				log.println(i+"\t"+aa+"\t"+uu+"\t"+xx+"\t"+bb+"\t"+fuu+"\t"+fxx);
			}
			if(bb-aa<gstol) {
				break;
			}
		}
		if(bb-aa>gstol) {
			log.println("Warning: Maximation did not converge within tolerance "+gstol+" in "+maxits+" iterations.");
		}
	}
	
	private void calcCI() {
		
		cthresh=fmle-1.92;
		double ul=getUpperLimit();
		double ll=getLowerLimit();
		CI[0]=ll;
		CI[1]=ul;
				
	}
	
	private double calcLL(double alpha) {
	
		double loglik=0;
		double loglik1=0;
		//log.println("*");
		for(SnpPair pr: paircts.keySet()) {
			HapCounter hapcts=paircts.get(pr);
			if(hapcts.totalCount()>0) {
				if (precise) {
					loglik+=hapcts.getLogLikMoreLogs(alpha);
				} else {
					loglik+=hapcts.getLogLik(alpha);
                }
            }
        }
		return(loglik);
	}
	
	private double[] bracketCI() {
	
		double minlesskey=0.5;
		double minlessval=ff_onehalf;
		double maxgreaterkey=mle;
		double maxgreaterval=fmle;
		
		for(Double key:ff.keySet()) {
			if(key > mle) {
				double val=ff.get(key);
				if(val<cthresh && val > minlessval) {
					minlesskey=key;
					minlessval=val;
				} else if (val > cthresh && val < maxgreaterval) {
					maxgreaterkey=key;
					maxgreaterval=val;
				}
			}
		}

		double maxlesskey=0;
		double maxlessval=ff_0;
        double mingreaterkey=mle;
        double mingreaterval=fmle;
        for(Double key:ff.keySet()) {
            if(key < mle) {
                double val=ff.get(key);
                if(val<cthresh && val > maxlessval) {
                    maxlesskey=key;
                    maxlessval=val;
                } else if (val > cthresh && val < mingreaterval) {
                    mingreaterkey=key;
                    mingreaterval=val;
                }
            }
        }
		if(debug) {
			log.println("conf_threshold="+cthresh+", lower bracket=("+maxlesskey+", "+mingreaterkey+"), ("+maxlessval+", "+mingreaterval+")");	
			log.println("conf_threshold="+cthresh+", upper bracket=("+minlesskey+", "+maxgreaterkey+"), ("+minlessval+", "+maxgreaterval+")");	
		}

		double[] ret={maxlesskey, mingreaterkey, minlesskey, maxgreaterkey};
		return ret;
	
	}
	
	private double getUpperLimit() {
	
		double minlesskey=0.5;
		double minlessval=ff_onehalf;
		double maxgreaterkey=mle;
		double maxgreaterval=fmle;
		for(Double key:ff.keySet()) {
			if(key > mle) {
				double val=ff.get(key);
				if(val<cthresh && val > minlessval) {
					minlesskey=key;
					minlessval=val;
				} else if (val > cthresh && val < maxgreaterval) {
					maxgreaterkey=key;
					maxgreaterval=val;
				}
			}
			
		}
		if(debug) {
			log.println("conf_threshold="+cthresh+", upper bracket=("+minlesskey+", "+maxgreaterkey+"), ("+minlessval+", "+maxgreaterval+")");
		}
		return bisect(maxgreaterkey, minlesskey, maxgreaterval, minlessval);
	
	}
	
	private double getLowerLimit() {
	
		double maxlesskey=0;
		double maxlessval=ff_0;
        double mingreaterkey=mle;
        double mingreaterval=fmle;
        for(Double key:ff.keySet()) {
            if(key < mle) {
                double val=ff.get(key);
                if(val<cthresh && val > maxlessval) {
                    maxlesskey=key;
                    maxlessval=val;
                } else if (val > cthresh && val < mingreaterval) {
                    mingreaterkey=key;
                    mingreaterval=val;
                }
            }
        }
		if(debug) {
			log.println("conf_threshold="+cthresh+", lower bracket=("+maxlesskey+", "+mingreaterkey+"), ("+maxlessval+", "+mingreaterval+")");	
		}
        return bisect(maxlesskey, mingreaterkey, maxlessval, mingreaterval);

	}
	
	private double bisect(double aa, double bb, double faa, double fbb) {
	
		double retval=-99;
		double xx=0.5*(aa+bb);
		
		if(bb-aa<gstol) {
			retval=xx;
		} else {
			double fxx=calcLL(xx);
			if(debug) {
				log.println("0\t"+aa+"\t"+xx+"\t"+bb+"\t"+fxx);
			}
			for(int i=1; i<maxits; i++) {
				double prod=(fxx-cthresh)*(fbb-cthresh);
				if(prod>0) {  // so fxx is on same side as fbb
					bb=xx;
					fbb=fxx;
				} else {
					aa=xx;
					faa=fxx;
				}
				xx=0.5*(aa+bb);
				fxx=calcLL(xx);
				if(debug) {
					log.println(i+"\t"+aa+"\t"+xx+"\t"+bb+"\t"+fxx);
				}
				retval=xx;
				if(bb-aa<gstol) {
					break;
				}
			}
			if(bb-aa >= gstol) {
				log.println("Root finding did not converge in within tolerance "+gstol+" in "+maxits+" iterations.");
			}
		}
		return retval;
	
	}
	
}

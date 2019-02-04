import org.apache.commons.math3.distribution.LogNormalDistribution;

import java.util.ArrayList;
import java.util.Random;

public class Microhabitat {

    //stuff for lognormal distribution of MICs
    private double mu = Math.log(7.92016113), sigma = 0.10018864;
    private LogNormalDistribution MIC_distribution = new LogNormalDistribution(mu, sigma);

    Random rand = new Random();
    private int K; //karrying kapacity
    private double c; //concentration of antimicrobial
    private double b; //migration rate for bacteria in this microhabitat
    ArrayList<Double> population; //contains a list of doubles representing the MICs of the bacteria present

    private boolean surface = false; //boolean detailing whether this is the microhabitat at the surface of the ship or not (maybe surface of biofilm too)
    private boolean biofilm_region = false; //boolean detailing whther this microhabitat is the well-mixed ocean environment or the biofilm region
    private double threshold_stickiness = 0.15; //fraction occupied required for biofilm formation to begin



    public Microhabitat(int K, double c){
        this.K = K;
        this.c = c;
        this.b = 0.1;
        this.population = new ArrayList<>(K);
    }


    public int getN(){return population.size();}
    public double getC(){return c;}

    public ArrayList<Double> getPopulation(){return population;}

    public void setSurface(boolean surface){this.surface = surface;}
    public boolean getBiofilm_region(){return this.biofilm_region;}
    public void setBiofilm_region(boolean biofilm_region){this.biofilm_region = biofilm_region;}
    public double getThreshold_stickiness(){return threshold_stickiness;}

    public double fractionFull(){
        return (double)getN()/(double)K;
    }


    public double stickiness(){
        //this returns a paramteter between 0 and 1 which gives the stickiness of the biofilm
        //the closer to 1, the stickier the biofilm, the less easy it is for the bacteria to move
        if(surface) return 1.;
        else if(getN() < threshold_stickiness*K) return 0.; //this is to prevent poisson errors

        else{
            double alpha = Math.log(2.)/((1. - threshold_stickiness)*K);
            return Math.min(Math.exp(alpha*(getN() - threshold_stickiness*K)) - 1., 1.);
        }
    }


    public double migrate_rate(){
        return b*(1. - stickiness());
    }


    public double beta(int index){
        return population.get(index);
    }

    public double phi_c(int index){
        double cB = c/beta(index);
        return 1. - 6*cB*cB/(5. + cB*cB);
    }

    public double replicationOrDeathRate(int index){
        //returns either +ve phi(c)*(1-N/k) [replication] or -ve phi(c) [death]
        //^TODO handle the double negative stuff here
        return  (phi_c(index) > 0.) ? phi_c(index)*(1. - getN()/K) : phi_c(index);
    }

    public double deathRate(int index){
        return phi_c(index);
    }




    public double getAvgGenotype(){
        if(getN() == 0) return 0.;
        else{
            double sum = 0.;
            for(Double geno : population){
                sum += geno;
            }
            return sum/(double)getN();
        }
    }


    public double getStDevOfGenotype(){
        if(getN() == 0) return 0.;
        else {
            double mean = getAvgGenotype();
            double sumSq = 0.;

            for(Double geno : population){
                sumSq += (geno-mean)*(geno-mean);
            }

            return Math.sqrt(sumSq/(getN()-1));
        }



    }



    public void randomlyPopulate(){

        for(int i = 0; i < K; i++){
            population.add(MIC_distribution.sample());
        }
    }

    public void randomlyPopulate(int randPopSize){
        for(int i = 0; i < randPopSize; i++){
            population.add(MIC_distribution.sample());
        }
    }


    public void addARandomBacterium(){
        population.add(MIC_distribution.sample());
    }

    public void replicateABacterium(int index){
        population.add(population.get(index));
    }

    public void replicateABacterium_x_N(int index, int nReps){
        for(int i = 0; i < nReps; i++){
            population.add(population.get(index));
        }
    }

    public void addABacterium(double MIC){population.add(MIC);}

    public void removeABacterium(int index){
        population.remove(index);
    }





}

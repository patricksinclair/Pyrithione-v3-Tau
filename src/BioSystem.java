import com.sun.xml.internal.org.jvnet.mimepull.MIMEConfig;
import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.Random;

public class BioSystem {

    Random rand = new Random();

    private int L, K;
    private double alpha, c_max;
    private Microhabitat[] microhabitats;
    private int initialMaxRandPop = 500; // I'll fill the first few microhabitats with a random number (0-maxrandpop) of random bacteria
    private int initialPopZone = 200;
    private double timeElapsed;
    private double tau;

    //counters to keep track of the number of events that happen

    private int deathCounter, replicationCounter, immigrationCounter,
            forcedImmigrationCounter, migrationCounter, nothingCounter;


    public BioSystem(int L, int K, double alpha, double c_max, double tau){
        this.L = L;
        this.K = K;
        this.alpha = alpha;
        this.tau = tau;
        this.timeElapsed = 0.;

        this.microhabitats = new Microhabitat[L];

        for(int i = 0; i < L; i++) {
            microhabitats[i] = new Microhabitat(K, BioSystem.getCValWithOffset(i, c_max, alpha, L));
        }

        microhabitats[L-1].setSurface(true);
        microhabitats[L-1].setBiofilm_region(true);
        microhabitats[L-1].randomlyPopulate(100);

        deathCounter = 0; replicationCounter = 0; immigrationCounter = 0;
        forcedImmigrationCounter = 0; migrationCounter = 0; nothingCounter = 0;

    }

    public double getTimeElapsed(){return timeElapsed;}
    public int getDeathCounter(){return deathCounter;}
    public int getReplicationCounter(){return replicationCounter;}
    public int getImmigrationCounter(){return immigrationCounter;}
    public int getForcedImmigrationCounter(){return forcedImmigrationCounter;}
    public int getMigrationCounter(){return migrationCounter;}
    public int getNothingCounter(){return nothingCounter;}


    public int getTotalN(){
        int runningTotal = 0;

        for(Microhabitat m : microhabitats){
            runningTotal += m.getN();
        }
        return  runningTotal;
    }


    public boolean everythingIsDead(){
        //returns true if the population of the system is 0

        return getTotalN() == 0;
    }

    public int getBiofilmEdge(){
        //finds the microhabitat furthest from the surface which is part of the biofilm
        int edgeIndex = 0;
        for(int i = 0; i < L; i++){
            if(microhabitats[i].getBiofilm_region()){
                edgeIndex = i;
                break;
            }
        }
        return edgeIndex;
    }

    public int getBiofilmSize(){
        // returns the no. of microhabitats which are in the biofilm block
        return L - getBiofilmEdge();
    }

    public int[] getPopulationDistribution(){
        int[] popSizes = new int[L];

        for(int i = 0; i < L; i++){
            popSizes[i] = microhabitats[i].getN();
        }
        return popSizes;
    }

    public double[] getAvgGenotypeDistribution(){
        double[] avgGenos = new double[L];
        for(int i = 0; i < L; i++){
            avgGenos[i] = microhabitats[i].getAvgGenotype();
        }
        return avgGenos;
    }



    public int[] getRandIndexes(int randBacteriaNumber){
        int totalCounter = 0;
        int microhab_index = 0;
        int bacteria_index = 0;

        int N_tot = getTotalN();
        // if we pick the bacteria 'outside the system', then we'll see if it gets immigrated
        if(randBacteriaNumber == N_tot) {
            return new int[]{-1, -1};

        } else {

            forloop:
            for(int i = 0; i < L; i++) {
                if(totalCounter + microhabitats[i].getN() <= randBacteriaNumber) {
                    totalCounter += microhabitats[i].getN();
                    continue forloop;

                } else {
                    microhab_index = i;
                    bacteria_index = randBacteriaNumber - totalCounter;
                    break forloop;
                }
            }
            return new int[]{microhab_index, bacteria_index};
        }
    }


    public void migrate(int mh_index, int bac_index){
        //Here this should only move bacteria within the biofilm region
        //takes the arguments of the microhabitat and bacteria indexes, so it performs
        //the necessary action on the desired one
        int biof_edge = getBiofilmEdge();
        double rand_bac = microhabitats[mh_index].getPopulation().get(bac_index);

        if(mh_index == L-1){
            microhabitats[mh_index].removeABacterium(bac_index);
            microhabitats[mh_index-1].addABacterium(rand_bac);

        }else if(mh_index == biof_edge){
            //if the bacteria is at the edge of the biofilm, there's a chance it detached, depending on the stickiness
            //ADD IN LATER
            microhabitats[mh_index].removeABacterium(bac_index);
            microhabitats[mh_index+1].addABacterium(rand_bac);

        }else{
            if(rand.nextBoolean()){
                microhabitats[mh_index].removeABacterium(bac_index);
                microhabitats[mh_index+1].addABacterium(rand_bac);
            }else{
                microhabitats[mh_index].removeABacterium(bac_index);
                microhabitats[mh_index-1].addABacterium(rand_bac);
            }
        }
    }


    public void replicate(int mh_index, int bac_index){
        microhabitats[mh_index].replicateABacterium(bac_index);
    }

    public void die(int mh_index, int bac_index){
        microhabitats[mh_index].removeABacterium(bac_index);
    }

    public void immigrate(int n_immigrants){
        for(int i = 0; i < n_immigrants; i++){
            microhabitats[getBiofilmEdge()-1].addARandomBacterium();
        }

    }


    public void migrate_v2(Microhabitat[] microhabitats, int mh_index, int nMigr, int originalPopSize){
        //performs a migration event on the designated microhabitat a designated no. of times
        //sends a random bacteria in a random direction
        int biof_edge = getBiofilmEdge();

        for(int i = 0; i < nMigr; i++) {

            int bac_index = rand.nextInt(originalPopSize);
            double rand_bac = microhabitats[mh_index].getPopulation().get(bac_index);

            if(mh_index == L - 1) {
                microhabitats[mh_index].removeABacterium(bac_index);
                microhabitats[mh_index - 1].addABacterium(rand_bac);

            } else if(mh_index == biof_edge) {
                //if the bacteria is at the edge of the biofilm, there's a chance it detached, depending on the stickiness
                //ADD IN LATER
                microhabitats[mh_index].removeABacterium(bac_index);
                microhabitats[mh_index + 1].addABacterium(rand_bac);

            } else {
                if(rand.nextBoolean()) {
                    microhabitats[mh_index].removeABacterium(bac_index);
                    microhabitats[mh_index + 1].addABacterium(rand_bac);
                } else {
                    microhabitats[mh_index].removeABacterium(bac_index);
                    microhabitats[mh_index - 1].addABacterium(rand_bac);
                }
            }
        }
    }


    public void replicate_v2(int mh_index, int bac_index, int nReps){

    }



    public void updateBiofilmSize(){
        //TODO update this so you don't get regions of non-biofilm in between biofilm regions
        for(Microhabitat m : microhabitats){
            if(m.fractionFull() >= m.getThreshold_stickiness()) m.setBiofilm_region(true);
        }
    }


    public void performAction(){
        // tau leaping implementation
        // go along each bacteria/microhabitat and work out the rates
        // then use these values as the means in a poisson distribution to work out how many times these events happen
        // array of ints corresponding to the number of events
        // for each microhabitat, have an array of ints corresponding to no. of migrations in/out
        // ?? get total no. of migrations, then randomly assign which direction they go ??
        // for each bacteria, an array of arrays for the no. of replications in each microhabitat for each species
        // ?? same with death ??
        // ??should the immigrants be added at the start or end of the event??
        // clear all extraneous bacteria from the system if they get detached maybe
        // ?? make new system, work out updates then replace the old one, or just iterate along ??
        // need to get replication rates before any replication is done, as otherwise the early replics will affect later ones
        // should we not allow bacteria that migrate into a microhab to be migrated again on the same iteration

        Microhabitat[] updated_microhabs = microhabitats.clone();

        double immigrationRate  = 100.;
        int nImmigrants = new PoissonDistribution(immigrationRate).sample();

        int biofilmSize = getBiofilmSize();
        int[][] replicationAllocations = new int[biofilmSize][];
        int[][] deathAllocations = new int[biofilmSize][];
        int[] migrationAllocations = new int[biofilmSize];



        int mh_counter = 0;
        for(int mh_index = getBiofilmEdge(); mh_index < L; mh_index++){

            int mh_pop = microhabitats[mh_index].getN();

            int n_migrations = new PoissonDistribution(microhabitats[mh_index].stickiness()).sample();

            int[] n_replications = new int[mh_pop];
            int[] n_deaths = new int[mh_pop];

            for(int bac_index = 0; bac_index < mh_pop; bac_index++){
                //works out the no. of replications each bacteria undergoes
                double gOrDRate_i = microhabitats[mh_index].replicationOrDeathRate(bac_index);

                if(gOrDRate_i >= 0){
                    n_replications[bac_index] = new PoissonDistribution(gOrDRate_i).sample();
                    n_deaths[bac_index] = 0;
                }
                else {
                    n_replications[bac_index] = 0;
                    n_deaths[bac_index] = new PoissonDistribution(Math.abs(gOrDRate_i)).sample();
                }
            }

            replicationAllocations[mh_counter] = n_replications;
            deathAllocations[mh_counter] = n_replications;
            migrationAllocations[mh_counter] = n_migrations;

            mh_counter++;
        }


        int mh_counter2 = 0;
        for(int mh_index = getBiofilmEdge(); mh_index < L; mh_index++){

            System.out.println("mh_index = "+mh_index);

            int originalPopSize = microhabitats[mh_index].getN(); //the original length of the arraylist before anything happens

            //iterate backwards so we can also remove bacteria
            for(int bac_index = originalPopSize-1; bac_index >= 0; bac_index--){
                System.out.println("bac_index: "+bac_index+"\toriginalPopSize: "+originalPopSize);
                updated_microhabs[mh_index].replicateABacterium_x_N(bac_index, replicationAllocations[mh_counter2][bac_index]);

                if(deathAllocations[mh_counter2][bac_index] > 0) updated_microhabs[mh_index].removeABacterium(bac_index);
            }

            //this needs to be done on updated_microhabitats
            System.out.println("original pop size: "+originalPopSize);
            System.out.println("nMigrations: "+ migrationAllocations[mh_counter2]);
            migrate_v2(updated_microhabs, mh_index, migrationAllocations[mh_counter2], originalPopSize);

        }

        //replace old system with new
        microhabitats = updated_microhabs;
        immigrate(nImmigrants);
        updateBiofilmSize();
        timeElapsed += tau;
    }



    public static void tester(){

        double duration = 200;
        int L = 500;
        int K = 500;
        double c_max = 10.;
        double alpha = 0.01;
        double tau = 0.01;

        int nTimeMeasurements = 20;
        double interval = duration/(double)nTimeMeasurements;

        BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);

        boolean alreadyRecorded = false;
        while(bs.timeElapsed <= duration){

           // if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.1*interval) && !alreadyRecorded){

                System.out.println("\ntime elapsed: "+String.valueOf(bs.getTimeElapsed())+"\ttotal N: "+String.valueOf(bs.getTotalN()));

                //System.out.println(Arrays.toString(bs.getNutrientsArray()));
                alreadyRecorded = true;

            //}
            //if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;
            bs.performAction();
        }
    }









    public static void getNumberOfEvents(double alpha, double tau){

        int K = 500, L = 500;
        int nReps = 10;
        int duration = 10000;
        int nCounters = 6;
        double c_max = 10.;
        double interval = duration/20.;
        Random rando = new Random();

        String filename = "pyrithione-randAlgorithmEvents";

        int[][] allEventCounts = new int[nReps][];

        String[] counterHeaders = {"migration", "immigration", "forced immigr", "replication", "death", "nothing"};


        for(int r = 0; r < nReps; r++){

            BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);

            while(bs.getTimeElapsed() <= duration){

                bs.performAction();
                //System.out.println("test");
                double timeElapsed = bs.getTimeElapsed();

                //if(timeElapsed%interval >= 0. && timeElapsed%interval <= 0.01){
                if(rando.nextInt(100)==50){
                    System.out.println("rep: "+String.valueOf(r)+"\ttime: "+String.valueOf(timeElapsed)+"\tN: "+String.valueOf(bs.getTotalN()));
                }
            }

            int[] eventCounts = {bs.getMigrationCounter(), bs.getImmigrationCounter(), bs.getForcedImmigrationCounter(),
                    bs.getReplicationCounter(), bs.getDeathCounter(), bs.getNothingCounter()};

            allEventCounts[r] = eventCounts;
        }

        double[] avgResults = Toolbox.averagedResults(allEventCounts);

        Toolbox.writeSimpleArrayAndHeadersToFile(filename, counterHeaders, avgResults);


    }














    public static double getCValWithOffset(int index, double maxC, double alpha, int L){
        //this calculates i* for the gradient profile offset, moves so the final concn is maxC, and it decreases with 1/e
        //or something like that
        //then calculates the corresponding concentration in that microhabitat

        double offset =  (L-1.) - Math.log(maxC+1.)/alpha;

        return (index >= offset) ? Math.exp(alpha*(index - offset)) - 1. : 0.;
    }
}

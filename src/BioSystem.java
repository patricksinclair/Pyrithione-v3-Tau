import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.Random;

public class BioSystem {

    Random rand = new Random();

    private int L, K;
    private double alpha, c_max;
    private Microhabitat[] microhabitats;
    private double timeElapsed;
    private double tau;

    //counters to keep track of the number of events that happen
    private int doubleDeathCounter;


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
        microhabitats[L-1].randomlyPopulate(25);

        doubleDeathCounter = 0;
    }

    public double getTimeElapsed(){return timeElapsed;}
    public int getDoubleDeathCounter(){return doubleDeathCounter;}


    public int getTotalN(){
        int runningTotal = 0;

        for(Microhabitat m : microhabitats){
            runningTotal += m.getN();
        }
        return  runningTotal;
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

    public double[] getPopulationDistribution(){
        double[] popSizes = new double[L];

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

    public double[] getStDevOfGenotypeDistribution(){
        double[] genoStDevs = new double[L];
        for(int i = 0; i < L; i++){
            genoStDevs[i] = microhabitats[i].getStDevOfGenotype();
        }
        return genoStDevs;
    }

    public int[][] getXountersArray(String[] headers){
        int[][] countersArray = new int[headers.length][L];
        for(int i = 0; i < L; i++){
            countersArray[0][i] = microhabitats[i].getImmigrationCounter();
            countersArray[1][i] = microhabitats[i].getMigrationInCounter();
            countersArray[2][i] = microhabitats[i].getMigrationOutCounter();
            countersArray[3][i] = microhabitats[i].getReplicationCounter();
            countersArray[4][i] = microhabitats[i].getDeathCounter();
        }
        return countersArray;
    }


    public void replicate(int mh_index, int bac_index){
        microhabitats[mh_index].replicateABacterium(bac_index);
    }

    public void die(int mh_index, int bac_index){
        microhabitats[mh_index].removeABacterium(bac_index);
    }

    public void immigrate(int n_immigrants, int mh_index){
        for(int i = 0; i < n_immigrants; i++){
            microhabitats[mh_index].addARandomBacterium();
        }

    }


    public void migrate_v2(Microhabitat[] microhabitats, int mh_index, int nMigr, int originalPopSize){
        //performs a migration event on the designated microhabitat a designated no. of times
        //sends a random bacteria in a random direction
        int biof_edge = getBiofilmEdge();

        //TODO this stops migrations being done on empty microhabs
        if(originalPopSize > 0) {
            for(int i = 0; i < nMigr; i++) {

                int bac_index = rand.nextInt(originalPopSize);
                double rand_bac = microhabitats[mh_index].getPopulation().get(bac_index);

                if(mh_index == L - 1) {
                    microhabitats[mh_index].removeABacterium(bac_index);
                    microhabitats[mh_index - 1].addABacterium(rand_bac);
                    microhabitats[mh_index-1].updateMigrationInCounter(1);

                } else if(mh_index == biof_edge) {
                    //if the bacteria is at the edge of the biofilm, there's a chance it detached, depending on the stickiness
                    //ADD IN LATER
                    microhabitats[mh_index].removeABacterium(bac_index);
                    microhabitats[mh_index + 1].addABacterium(rand_bac);
                    microhabitats[mh_index+1].updateMigrationInCounter(1);

                } else {
                    if(rand.nextBoolean()) {
                        microhabitats[mh_index].removeABacterium(bac_index);
                        microhabitats[mh_index + 1].addABacterium(rand_bac);
                        microhabitats[mh_index+1].updateMigrationInCounter(1);
                    } else {
                        microhabitats[mh_index].removeABacterium(bac_index);
                        microhabitats[mh_index - 1].addABacterium(rand_bac);
                        microhabitats[mh_index-1].updateMigrationInCounter(1);
                    }
                }
            }
        }
    }



    public void updateBiofilmSize(){
        //TODO update this so you don't get regions of non-biofilm in between biofilm regions
        //actually that can't happen, as we never set it to false in this loop
        //IMMOD1 changed the for loop so that the biofilm edge is the empty microhab at the end
        /*for(Microhabitat m : microhabitats){
            if(m.fractionFull() >= m.getThreshold_stickiness()) m.setBiofilm_region(true);
        }*/
        forloop:
        for(int i=L-1; i > 0; i--){
            if(microhabitats[i].fractionFull() >= microhabitats[i].getThreshold_stickiness()){
                microhabitats[i].setBiofilm_region(true);

                if(microhabitats[i-1].fractionFull()<microhabitats[i].getThreshold_stickiness()){
                    microhabitats[i-1].setBiofilm_region(true);
                    break forloop;
                }
            }
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
        int nImmigrants = new PoissonDistribution(immigrationRate*tau).sample();

        int biofilmSize = getBiofilmSize();
        int[][] replicationAllocations = new int[biofilmSize][];
        int[][] deathAllocations = new int[biofilmSize][];
        int[] migrationAllocations = new int[biofilmSize];


        int mh_counter = 0;
        for(int mh_index = getBiofilmEdge(); mh_index < L; mh_index++){

            int mh_pop = microhabitats[mh_index].getN();

            double migration_rate = microhabitats[mh_index].migrate_rate();
            int n_migrations = (migration_rate > 0) ? new PoissonDistribution(migration_rate*tau).sample() : 0;

            int[] n_replications = new int[mh_pop];
            int[] n_deaths = new int[mh_pop];

            for(int bac_index = 0; bac_index < mh_pop; bac_index++){
                //works out the no. of replications each bacteria undergoes
                double gOrDRate_i = microhabitats[mh_index].replicationOrDeathRate(bac_index);

                if(gOrDRate_i == 0){
                    //this stops the poisson distribution getting funky with mean=0
                    n_replications[bac_index] = 0;
                    n_deaths[bac_index] = 0;

                }else if(gOrDRate_i > 0 ){
                    n_replications[bac_index] = new PoissonDistribution(gOrDRate_i*tau).sample();
                    n_deaths[bac_index] = 0;
                }
                else {
                    n_replications[bac_index] = 0;
                    n_deaths[bac_index] = new PoissonDistribution(Math.abs(gOrDRate_i)*tau).sample();
                }
            }

            replicationAllocations[mh_counter] = n_replications;
            deathAllocations[mh_counter] = n_deaths;
            migrationAllocations[mh_counter] = n_migrations;

            mh_counter++;
        }


        int mh_counter2 = 0;
        for(int mh_index = getBiofilmEdge(); mh_index < L; mh_index++){

            int originalPopSize = microhabitats[mh_index].getN(); //the original length of the arraylist before anything happens

            //iterate backwards so we can also remove bacteria
            for(int bac_index = originalPopSize-1; bac_index >= 0; bac_index--){

                updated_microhabs[mh_index].replicateABacterium_x_N(bac_index, replicationAllocations[mh_counter2][bac_index]);
                updated_microhabs[mh_index].updateReplicationCounter(replicationAllocations[mh_counter2][bac_index]);

                if(deathAllocations[mh_counter2][bac_index] > 0) updated_microhabs[mh_index].removeABacterium(bac_index);
                if(deathAllocations[mh_counter2][bac_index] > 1) doubleDeathCounter++;
                updated_microhabs[mh_index].updateDeathCounter(deathAllocations[mh_counter2][bac_index]);

            }

            mh_counter2++;
        }

        int mh_counter_migration = 0;
        for(int mh_index = getBiofilmEdge(); mh_index < L; mh_index++){

            int originalPopSize = microhabitats[mh_index].getN();

            migrate_v2(updated_microhabs, mh_index, migrationAllocations[mh_counter_migration], originalPopSize);
            updated_microhabs[mh_index].updateMigrationOutCounter(migrationAllocations[mh_counter_migration]);
            mh_counter_migration++;
        }

        updated_microhabs[getBiofilmEdge()].updateImmigrationCounter(nImmigrants);
        //replace old system with new
        microhabitats = updated_microhabs;
        //IMMOD1 removed the -1 from immigration index
        immigrate(nImmigrants, getBiofilmEdge());
        updateBiofilmSize();
        timeElapsed += tau;
    }



    public static void getPopDistbInfo(){
        //method to get info on population distbs
        //get popsize over time
        //pop distb over time
        //biofilm edge over time
        //avg genotype distb over time

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nReps = 10, nMeasurements = 20;
        double duration = 200., interval = duration/nMeasurements;

        String popSizeFilename = "pyrithione-testing-pop_size-t="+String.valueOf(duration);
        String popDistbFilename = "pyrithione-testing-pop_distb-t="+String.valueOf(duration);
        String biofilmEdgeFilename = "pyrithione-testing-biofilm_edge-t="+String.valueOf(duration);
        String avgGenotypeDistbFilename = "pyrithione-testing-avgGenoDistb-t="+String.valueOf(duration);
        String genoStDevDistbFilename = "pyrithione-testing-genoStDevDistb-t="+String.valueOf(duration);
        String counterDistbsFilename = "pyrithione-testing-counterDistb-t="+String.valueOf(duration);

        String[] counterHeaders = {"immigration", "migrationIn", "migrationOut", "replication", "death"};
        int nCounters = counterHeaders.length;

        double[][] allPopSizes = new double[nReps][];
        double[][][] allPopDistbs = new double[nReps][][];
        double[][] allBiofilmEdges = new double[nReps][];
        double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        double[][][] allGenoStDevs = new double[nReps][][];
        double[][] allCounters = new double[nReps][];


        for(int r = 0; r < nReps; r++){

            BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);

            boolean alreadyRecorded = false;
            int timerCounter = 0;

            double[] popSizes = new double[nMeasurements+1];
            double[][] popDistbs = new double[nMeasurements+1][];
            double[] biofilmEdges = new double[nMeasurements+1];
            double[][] avgGenotypeDistbs = new double[nMeasurements+1][];
            double[][] genoStDevs = new double[nMeasurements+1][];


            while(bs.timeElapsed <= duration+0.02*interval){

                if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.1*interval) && !alreadyRecorded){

                    String output = String.format("rep: %d \ttime elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d \tdouble deaths: %d",
                            r, bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getDoubleDeathCounter());

                    System.out.println(output);

                    popSizes[timerCounter] = bs.getTotalN();
                    popDistbs[timerCounter] = bs.getPopulationDistribution();
                    biofilmEdges[timerCounter] = bs.getBiofilmEdge();
                    avgGenotypeDistbs[timerCounter] = bs.getAvgGenotypeDistribution();
                    genoStDevs[timerCounter] = bs.getStDevOfGenotypeDistribution();

                    alreadyRecorded = true;
                    timerCounter++;
                }
                if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

                bs.performAction();
            }

            allPopSizes[r] = popSizes;
            allPopDistbs[r] = popDistbs;
            allBiofilmEdges[r] = biofilmEdges;
            allAvgGenotypeDistbs[r] = avgGenotypeDistbs;
            allGenoStDevs[r] = genoStDevs;
        }

        double[] processedPopSizes = Toolbox.averagedResults(allPopSizes);
        double[][] processedPopDistbs = Toolbox.averagedResults(allPopDistbs);
        double[] processedBiofilmEdges = Toolbox.averagedResults(allBiofilmEdges);
        double[][] processedAvgGenotypeDistbs = Toolbox.averagedResults(allAvgGenotypeDistbs);
        double[][] processedGenoStDevs = Toolbox.averagedResults(allGenoStDevs);

        Toolbox.writeAveragedArrayToFile(popSizeFilename, processedPopSizes);
        Toolbox.writeAveragedDistbsToFile(popDistbFilename, processedPopDistbs);
        Toolbox.writeAveragedArrayToFile(biofilmEdgeFilename, processedBiofilmEdges);
        Toolbox.writeAveragedDistbsToFile(avgGenotypeDistbFilename, processedAvgGenotypeDistbs);
        Toolbox.writeAveragedDistbsToFile(genoStDevDistbFilename, processedGenoStDevs);

        System.out.println("results written to file");
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

            String output = String.format("time elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d \tdouble deaths: %d",
                    bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getDoubleDeathCounter());

            System.out.println(output);

            //System.out.println(Arrays.toString(bs.getNutrientsArray()));
            alreadyRecorded = true;

            //}
            //if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;
            bs.performAction();
        }
    }








    public static double getCValWithOffset(int index, double maxC, double alpha, int L){
        //this calculates i* for the gradient profile offset, moves so the final concn is maxC, and it decreases with 1/e
        //or something like that
        //then calculates the corresponding concentration in that microhabitat

        double offset =  (L-1.) - Math.log(maxC+1.)/alpha;

        return (index >= offset) ? Math.exp(alpha*(index - offset)) - 1. : 0.;
    }
}

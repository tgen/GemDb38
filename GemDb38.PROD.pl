#!/usr/bin/perl -w
#/*************************************************************************
# *
# * TGen CONFIDENTIAL
# * __________________
# *
# *  [2010] - [2014] Translational Genomics Research Institute (TGen)
# *  All Rights Reserved.
# *
# * NOTICE:  All information contained herein is, and remains
# * the property of Translational Genomics Research Institute (TGen).
# * The intellectual and technical concepts contained herein are proprietary
# * to  TGen and may be covered by U.S. and Foreign Patents,
# * patents in process, and are protected by trade secret or copyright law.
# * Dissemination of this information, application or reproduction
# * is strictly forbidden unless prior written permission is obtained
# * from TGen.
# *
# * Major Contributor(s):
#    David Craig
#    Release 12/31/14
#
#  Dependencies
#     MongoDB,Time::localtime,File::stat, Scalar::Util, Storable, FindBin,
#/
#######   Determines which Configuration Library is needed ###########################


##### DCUNITS ######

our $runType = "PROD";
our $VERSION = '3.95';
##### DCUNITS ######


######################################################################################
###                     These should not need to be changed                       ####
###                      Please edit conf file                                    ####
######################################################################################
$| = 1;
our $GemStatus = 0;
use FindBin '$Bin';
use Time::localtime;
require "$Bin/lib/Maintain.$VERSION.pm";
require "$Bin/lib/InsertVar.$VERSION.pm";
require "$Bin/lib/Annotate.$VERSION.pm";
require "$Bin/lib/Heuristics.$VERSION.pm";
our $date          = localtime->year * 365 * 3600 * 24 + localtime->yday * 3600 * 24 + localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
our $start         = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
our $stdout        = *STDOUT;
our $stderr        = *STDERR;
our $diff          = $start-$start;
our $RunParameters = { 'runDir' => $Bin, 'version' => $VERSION, 'binDir' => $Bin };
chdir( $RunParameters->{'runDir'} );


######################################################################################    
###                      MAIN BLOCK                                               ####
######################################################################################

{  
##### HELP PRINTOUT######
if (join(@ARGV)=~/--h/) {
    print <<OUT
    _____________GemDb38 $VERSION, Copyright 2014 TGen, David W. Craig_____________________
    Requires "GemDb38.PROD.conf.pm" which contains run parameters (RunParameters)
    Example Commands:
         ./GemDb38.pl     
         ./GemDb38.pl threads=16   
         ./GemDb38.pl runType=SU2C    
         ./GemDb38.pl clean=1 
        
    Options:   
          clean=1       | Removes all partial or incomplete runs, and retries
          nolog=1       | No log files (print to STDOUT/STDERR      
          runType=SU2C  | Parameter file is GemDb38.SU2C.conf.txt, default runType is PROD
          CreateOnPath= | Full Path for Deletion of all contained VCFs, then adding VCFs
          threads=16    | Start with 16 threads
      
    Examples:
          ./GemDb38.pl                                  |  Rebuilds from scratch
          ./GemDb38.pl clean=1                          |  Cleans from crash
          ./GemDb38.pl AddOnly=su2c                     |  Adds specific files with 'su2c' in name  
          ./GemDb38.pl ForcePath=su2c runType=SU2C      |  Removes any entrys and reinserts with new runType    
          ./GemDb38.pl CreateOnPath=/reports/C017/0001/20140609/T3/A1STX/analysisResults/ 

    Files For Communication:
          PROD.RUNNING        |  At least one sample has been added.  If not present, database is purged.
          PROD.pid123.ADDING  |  Markers being added by process PID 123
          PROD.INACTIVE       |  No records are being added.    
          PROD.TOCLEAN        |  Flags that records have been added, and clean should start when nothing left to add   
          PROD.pid123.CLEANING|  Cleaning process is underway    
          PROD.INACTIVE       |  No records are being added.    
    -------END OF HELP-----


___________ STARTING GemDb38 $VERSION, Copyright 2014 TGen, David W. Craig ______
            #######-- Configuration File is GemDb38.$runType.conf.pm -#########
_______________________________________________________________________________   

OUT
      ;
}
###############################################################    
## Clean OLD files
###############################################################    
    unless (-e "$Bin/log") {system("mkdir $Bin/log")}
    unless (-e "$Bin/runDir") {system("mkdir $Bin/runDir")}
    system("find ./log/ -name \"GemDb38.$runType.*log\" -mtime +1 -exec cat {} >> ./log/$runType.old.log \\;");
    system("find ./log/ -name \"GemDb38.$runType.*log\" -mtime +1 -exec rm {} \\;");
    print "++---------------------------------------------++\n";
    print "+For help use --h.  Log printing to $RunParameters->{'runDir'}/log/GemDb38.$runType.$date.pid$$.log. To print log to STDOUT use -nolog=1.  \nTo suspend and push to background use control-z and 'bg;disown;'\n   To clean crashed runs, use ./GemDb.pl clean=1 or rm runDir/DCUNIT.pid*\n";
    print "\n++----------- GemDb38 Running -------------------------++\n";
    open my $LOG, ">", "$RunParameters->{'runDir'}/log/GemDb38.$runType.$date.pid$$.log" or die "Can't open $RunParameters->{'runDir'}/GemDb38.$date.log";
    $RunParameters->{'LOG'} = \$LOG;

###############################################################    
##### LOAD COMMANDLINE ARGS ######
###############################################################
  $RunParameters->{'threads'}=1;
  ARGS: foreach my $arg (@ARGV) {
        chomp($arg);
        if ( $arg =~ /(.*?)=(.*)/ ) {
            if ( length($1) > 0 && length($2) > 0 ) {
                chomp($2);
                $RunParameters->{$1} = $2;
                $RunParameters->{'StartupRunParameters'}->{$1} = $2;
            }
        }
    }
    if ( exists( $RunParameters->{runType} ) ) { $runType = $RunParameters->{runType}; }
    print "\n\n\n++------------------------------------------------++\n";
    print "++------------------------------------------------++\n";          
    print "++---GemDb38 $VERSION started $date and pid:$$---++\n";
    print "++------------------------------------------------++\n";              
    print "++------------------------------------------------++\n\n\n";    
    if ($RunParameters->{'threads'} > 1) {
        print "\t+Starting Threading\n";
        for ($j=1;$j<=$RunParameters->{'threads'};++$j) {
            sleep(1);
            $procCount=`ps -ef | grep GemDb | grep perl | wc -l`; chomp($procCount); 
            --$procCount;
            if ($procCount >= $RunParameters->{'threads'}) { 
                 print "\t$procCount Threads already started, $RunParameters->{'threads'} requested.\n";
                 goto END;
            }
            print "\t+Forking thread $j, $procCount threads running\n";
            system("nohup $Bin/GemDb38.$runType.pl >& /dev/null &");
        }
        goto END;
    }
    
    unless ( join( ":", @ARGV ) =~ /NoLog/i ) { *STDERR = $LOG; *STDOUT = $LOG; }
###############################################################
##### LOAD RunParamaters######
###############################################################
    print "\n\n #######---Configuration File is GemDb38.$runType.conf.pm--#########\n";
    system("rm -f $RunParameters->{'runDir'}/$runType.INACTIVE");
    if ( -e "$Bin/GemDb38.$runType.conf.pm" ) {
        require "$Bin/GemDb38.$runType.conf.pm";
        $started                = "$runType.RUNNING";
        $RunParameters          = Conf->loadDefaults($RunParameters);
        $RunParameters->{'LOG'} = \$LOG;
        foreach $arg (@ARGV) {
            chomp($arg);
            if ( $arg =~ /(.*?)=(.*)/ ) { $RunParameters->{$1} = $2 }
        }
        Maintain->mongoConnect($RunParameters);
    }
    else { die "!!!!Can't find GemDb38.$runType.conf.pm!!!\n"; }

    Maintain->printSystemDefaults($RunParameters);

###############################################################
##### STARTUP ######
###############################################################
    unless ( -e "$RunParameters->{'runDir'}/$started" ) {
        print "\n+++-------------- Replacing a DCUNIT.RUNNING file, did you mean to erase database.  Its still there-------------+++ \n";
        system("echo > $RunParameters->{'runDir'}/$started");
        $RunParameters->{'largeRun'} = int(1);
    #    Maintain->dropDb($RunParameters);
        $GemStatus = 1;
    }
###############################################################
##### CLEAN UP AFTER CRASH ######
###############################################################
    if ( exists( $RunParameters->{'purgePartials'} ) ) {
        print "\n+++-------------Purge Partial adds-------+++ \n";
        system("echo > $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
        $RunParameters->{'scrubEnabled'}=1;
        system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
        system("rm -f $RunParameters->{'runDir'}/$runType.TOCLEAN");
        system("rm -f $RunParameters->{'runDir'}/$runType.*.ADDING");    
        $GemStatus = 1;        
        Maintain->purgePartials($RunParameters);             
    }
    if ( exists( $RunParameters->{'clean'} ) ) {
        print "\n+++-------------Cleaning and Restarting Due to Crash-------+++ \n";
        system("echo > $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
        $RunParameters->{'scrubEnabled'}=1;
        system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
        system("rm -f $RunParameters->{'runDir'}/$runType.TOCLEAN");
        system("rm -f $RunParameters->{'runDir'}/$runType.*.ADDING");        
        Maintain->closeAllConnections($RunParameters);
        $GemStatus = 1;        
        Maintain->purgePartials($RunParameters);             
        #Maintain->cleanAndMerge($RunParameters);
        Maintain->calculateDbFreq($RunParameters);        
        Heuristics->indexCollections($RunParameters);
    }
    if ( glob("$RunParameters->{'runDir'}/$runType.pid*") ) { 
        $RunParameters->{'pidRunning'} = 1; 
    } else {
         $RunParameters->{'pidRunning'} = 0;
    }
    
###############################################################
#### INSERTING ######
###############################################################
    if ( -e "$RunParameters->{'runDir'}/$started" && !( glob("$RunParameters->{'runDir'}/$runType.pid*.CLEANING") ) ) {
        print "\n+++--------------GemDb38 Adding Files--------------+++ \n";
        system("echo > $RunParameters->{'runDir'}/$runType.pid$$.ADDING");
        Maintain->buildBufferConnections($RunParameters);
        if ($RunParameters->{'pidRunning'} ==0) {
           if (exists($RunParameters->{'sync'} ) ) {
               if ($RunParameters->{'sync'}==1) {
                  $RunParameters->{'clean'}=1;
               }
           }        
        }
        if ( InsertVar->Insert($RunParameters) > 0 ) {
            Annotate->annotate( $RunParameters);
            print("Setting GemStatus to 1 in Annotate\n");
	    $GemStatus = 1;
            system("echo > $RunParameters->{'runDir'}/$runType.TOCLEAN");
        }
        system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.ADDING");
        Maintain->closeBufferConnections($RunParameters);
    }

###############################################################
##### CLEANING MAINTAINANCE & FIRST RUN ######
###############################################################

    unless ( glob("$RunParameters->{'runDir'}/$runType.pid*.ADDING") || glob("$RunParameters->{'runDir'}/$runType.pid*.CLEANING") ) {
        if ( -e "$RunParameters->{'runDir'}/$runType.TOCLEAN" ) {
            print "\n+++--------------GemDb38 Clean, Merge, Indexing--------------+++ \n";
            system("echo > $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
            Maintain->calculateDbFreq($RunParameters);
            Heuristics->indexCollections($RunParameters);
            Heuristics->runInheritance($RunParameters,'germline');
            Heuristics->joinVCF($RunParameters);
            $status=1;
#            $status=Maintain->cleanAndMerge($RunParameters);                        
            $GemStatus = 1;
            system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
            if ($status==1) {
                system("rm -f $RunParameters->{'runDir'}/$runType.TOCLEAN");
            }
        }
        else {
            print "\t+Nothing to Clean\n";
        }
    }
    system("echo > $RunParameters->{'runDir'}/$runType.INACTIVE");


###############################################################
##### SPAWN IF NOT FINISHED ######
###############################################################
    $diff = sprintf( "%3.2f", ( localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start ) / 60 );
    if ( $diff < 0 ) { $diff = $diff + 1440 };
    print "\n+++-----GemDb38 Finished--($diff minutes)--------+++ \n";
    close $LOG;
}

#####  EXIT STRATEGY  #####
sleep(2);
END: if ( $GemStatus == 0 ) {
    system("rm  ./log/GemDb38.$runType.$date.pid$$.log");
    print "--Ending Gracefully--\n";
} elsif ( $diff > 5 ) {
    system("cat ./log/$runType.old.log ./log/GemDb38.$runType.$date.pid$$.log > ./log/$runType.old.log");
    system("rm  ./log/GemDb38.$runType.$date.pid$$.log");
    print $stdout "\t+Spawning new process\n";
    system("$Bin/GemDb38.$runType.pl&");
    print "--Ending This Process pid$$ Gracefully--\n";
}

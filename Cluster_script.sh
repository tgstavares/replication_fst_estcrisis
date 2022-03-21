#!/bin/bash

### Simple script to create clusters in AWS with mpi capabilities

# Cluster characteristics
NUM_INST=81
NAME_CLT="TiagoCluster"

sec_key="" # choose an appropriate security key
ami_omp="" # choose an appropriate ami
ami_mpi="" # choose an appropriate ami
ami_acc="" # choose an appropriate ami
int_omp="c5a.4xlarge"
int_mpi="c5.large" # t2.micro; t3.nano; t3a.micro; t2.medium; c5a.large; t1.micro; c5.large;; c6g.large (arm)
int_acc="g4dn.xlarge"

if [ -z $1 ];then
    ## ERROR MESSAGES
    echo
    echo "  ERROR: Input the first argument: start-spot, start-ondemand, or terminate"
    echo "  Usage: Cluster_script (start-spot | start-ondemand | reconfigure | terminate) [omp | acc | mpi]"
    echo
    exit    
elif [ "$1" = "start-spot" -o "$1" = "start-ondemand" -o "$1" = "reconfigure" -o "$1" = "terminate" ];then

    if [ "$1" = "start-spot" -o "$1" = "start-ondemand" ];then
	FILE=.Info_instances.txt
	if test -f "$FILE"; then
    	    TAG=$(awk 'FNR == 1 {print $6}' .Info_instances.txt)
    	    echo
    	    echo "  ERROR: Attempt to create a cluster with a $FILE; check for status; tag is $TAG"
    	    echo
	    exit
	else
	    if [ -z $2 ];then
		echo
		echo "  ERROR: Input the second argument: omp,acc, or mpi"
		echo "  Usage: Cluster_script (start-spot | start-ondemand | reconfigure | terminate) [omp | acc | mpi]"
		echo
		exit
	    elif [ "$2" = "omp" ];then
		AMI_SPEC="$ami_omp"
		TYP_INST="$int_omp"
		NUM_INST=1
		WTIME=60
	    elif [ "$2" = "mpi" ];then
		AMI_SPEC="$ami_mpi"
		TYP_INST="$int_mpi"
		WTIME=120
	    elif [ "$2" = "acc" ];then
		AMI_SPEC="$ami_acc"
		TYP_INST="$int_acc"
		NUM_INST=1
		WTIME=180
	    else
		echo
		echo "  ERROR: Input a valid second argument: omp,acc, or mpi"
		echo "  Usage: Cluster_script (start-spot | start-ondemand | reconfigure | terminate) [omp | acc | mpi]"
		echo
		exit
	    fi	    
	    
    	    start=$SECONDS

    	    ## Start cluster
    	    TAG_INST=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 6 | head -n 1)

    	    # Start script here
    	    echo
    	    echo "  Setting up $NAME_CLT cluster with $NUM_INST $TYP_INST instances..."
    	    echo

    	    # Launch instances - some of these things are specific like the ami, the key-name, and the security group -> change this at will	
    	    if [ "$1" = "start-ondemand" ];then
    		aws ec2 run-instances --image-id $AMI_SPEC --instance-type $TYP_INST --key-name mykey --security-group-ids $sec_key  --count $NUM_INST --tag-specifications ResourceType=instance,Tags=[\{Key=$NAME_CLT,Value=$TAG_INST\}] > /dev/null
    	    elif [ "$1" = "start-spot" ];then
    		aws ec2 run-instances --instance-market-options '{"MarketType":"spot"}' --image-id $AMI_SPEC --instance-type $TYP_INST --key-name mykey --security-group-ids $sec_key --count $NUM_INST --tag-specifications ResourceType=instance,Tags=[\{Key=$NAME_CLT,Value=$TAG_INST\}] > /dev/null
    	    else
    		echo
    		echo "  ERROR: Choose spot or ondemand pricing"
    		echo "  Usage: Cluster_script (start-spot | start-ondemand | reconfigure | terminate) [omp | acc | mpi]"
    		echo
    		exit
    	    fi	    

    	    i=1
    	    sp="/-\|"
    	    iter=0
    	    # This is mostly to wait until they are in a running state
    	    trap '[ -z $! ] || kill $!' SIGHUP SIGINT SIGQUIT SIGTERM
    	    aws ec2 wait instance-running --filters "Name=tag:$NAME_CLT,Values=$TAG_INST" &
    	    while [ -e /proc/$! ]; do
    		if [ $iter -ge 300 ]; then
    		    iter=1000
    		    break
    		fi
    		((iter++))
    		aws ec2 describe-instances --filters "Name=instance-state-name,Values=running" "Name=tag:$NAME_CLT,Values=$TAG_INST" \
    		    --output text --query "Reservations[].Instances[].[InstanceId,PublicIpAddress,PrivateIpAddress,State.Name,[Tags[?Key=='$NAME_CLT'].Value] [0][0]]" > temp0.txt
    		ncurinst=$(awk '($4 == "running") {count++ } END { print count }' temp0.txt)
    		[ -z $ncurinst ] && ncurinst=0
    		printf "\r  [${sp:i++%${#sp}:1}] Launching instances: $ncurinst of $NUM_INST"
		duration=$(( SECONDS - start ))
		if [ $ncurinst -eq 0 -a $duration -gt 60 ];then
    		    echo
    		    echo "  ERROR: cannot create the first instance ($ncurinst and $duration); check output log"
    		    echo "  Usage: Cluster_script (start-spot | start-ondemand | reconfigure | terminate) [omp | acc | mpi]"
    		    echo
    		    exit
		fi
    		sleep 1
    	    done

    	    # This is to check that all are running and get a temp info file
    	    aws ec2 describe-instances --filters "Name=instance-state-name,Values=running" "Name=tag:$NAME_CLT,Values=$TAG_INST" \
    		--output text --query "Reservations[].Instances[].[InstanceId,PublicIpAddress,PrivateIpAddress,State.Name,[Tags[?Key=='$NAME_CLT'].Value] [0][0]]" > temp0.txt
    	    ncurinst=$(awk '($4 == "running") {count++ } END { print count }' temp0.txt)
    	    [ -z $ncurinst ] && ncurinst=0
    	    if [ $iter -lt 999 ]; then
    		printf "\r  [t] Launching instances: $ncurinst of $NUM_INST"
    		echo; echo
    	    else
    		printf "\r  [f] Launching instances: $ncurinst of $NUM_INST --> TIMEOUT ERROR"
    		echo; echo
		exit
    	    fi

    	    # Now create instances info file in a human readable format and a hosts file with ips
    	    awk '{printf "node%03g   ",NR-1; print $0,NR}' temp0.txt > .Info_instances.txt
    	    awk '{print $4,$1}' .Info_instances.txt > .hosts_info.txt
    	    awk '{print $1}' .hosts_info.txt > .hosts_file.txt
    	    awk 'NR == 1 {print $4}' .Info_instances.txt > .ip_master.txt
    	    rm temp0.txt

    	    # Change instance names
    	    cat .Info_instances.txt | xargs -l -n7 bash -c 'aws ec2 create-tags --resources $1 --tag "Key=Name,Value=$0"'
    	    sleep 10
    	    duration=$(( SECONDS - start ))
    	    if [ $duration -lt $WTIME ];then
    		echo "  Wait..."
    		echo
    		sleep $WTIME
    	    fi
	fi
    fi    

    if [ "$1" = "start-spot" -o "$1" = "start-ondemand" -o "$1" = "reconfigure" ];then
    #if [ "$1" = "reconfigure" ];then

	FILE=.Info_instances.txt
	if test -f "$FILE"; then
     	    IPMASTER=$(awk 'NR == 1 {print $4}' .Info_instances.txt)
	    echo
	    echo "  Now setting up passwordless ssh"
	    echo

    	    echo "########################################"
    	    #cat .Info_instances.txt | xargs -l -n7 bash -c 'ssh -4 -i ~/.ssh/mykey2.rsa -o StrictHostKeyChecking=no -o GSSAPIAuthentication=no ubuntu@$3 "sudo cp /home/ubuntu/.ssh/authorized_keys /root/.ssh" && echo "p1 $0"'
	    parallel-ssh -i -h .hosts_file.txt -l ubuntu -x "-o StrictHostKeyChecking=no -i /home/tgst/.ssh/mykey2.rsa" 'sudo cp /home/ubuntu/.ssh/authorized_keys /root/.ssh'
	    sleep 1

	    #cat .Info_instances.txt | xargs -l -n7 bash -c 'scp -i ~/.ssh/mykey2.rsa ~/.ssh/mykey2.rsa ubuntu@$3:~/.ssh'
            #cat .Info_instances.txt | xargs -l -n7 bash -c 'scp -i ~/.ssh/mykey2.rsa -o StrictHostKeyChecking=no -o GSSAPIAuthentication=no ~/.ssh/mykey2.rsa root@$3:~/.ssh'
	    parallel-scp -h .hosts_file.txt -l root -x "-o StrictHostKeyChecking=no -i /home/tgst/.ssh/mykey2.rsa" /home/tgst/.ssh/mykey2.rsa /root/.ssh	    
	    sleep 1
	    
            #cat .Info_instances.txt | xargs -l -n7 bash -c 'scp -i ~/.ssh/mykey2.rsa -o GSSAPIAuthentication=no ./.hosts_info.txt root@$3:/etc/hosts'
    	    #cat .Info_instances.txt | xargs -l -n7 bash -c 'ssh -i ~/.ssh/mykey2.rsa root@$3 "sudo cat ./.hosts_info.txt >> /etc/hosts"'
	    #cat .Info_instances.txt | xargs -l -n7 bash -c 'ssh -i ~/.ssh/mykey2.rsa root@$3 "sudo cat ./.hosts_info.txt >> /etc/hosts" && echo "p2 $0"'
	    parallel-scp -h .hosts_file.txt -l root -x "-o StrictHostKeyChecking=no -i /home/tgst/.ssh/mykey2.rsa" .hosts_info.txt /etc/hosts
	    
    	    ssh -i ~/.ssh/mykey2.rsa root@$IPMASTER "echo 'eval \`ssh-agent\` &> /dev/null && ssh-add ~/.ssh/mykey2.rsa &> /dev/null' >> .profile"
    	    ssh -i ~/.ssh/mykey2.rsa root@$IPMASTER "eval \`ssh-agent\` && ssh-add ~/.ssh/mykey2.rsa && cat /etc/hosts | grep -i '.' | grep -v '#\|:' | awk '{ print \$2}' | xargs -I{} ssh-copy-id -o StrictHostKeyChecking=no {}"
	    scp -i ~/.ssh/mykey2.rsa .hosts_file.txt root@$IPMASTER:/home/ubuntu
    	    echo "########################################"
    	    echo

	    if [ $NUM_INST -gt 1 ]; then
		echo "  Now setting up a NFS folder"
		ssh -i ~/.ssh/mykey2.rsa root@$IPMASTER "sudo echo '/home *(async,no_root_squash,no_subtree_check,rw)' >> /etc/exports && exportfs -a"
		echo "  Mounting server"
		sleep 1
		tail -n+2 .Info_instances.txt | xargs -l -n7 bash -c 'ssh -i ~/.ssh/mykey2.rsa -o GSSAPIAuthentication=no root@$3 "sudo echo "node000:/home /home nfs rw,exec,noauto 0 0" >> /etc/fstab && sudo mount -t nfs node000:/home /home" && echo "p3 $0"'
		echo
	    fi
	    
	    echo "  Done! Access master using: ssh -i ~/.ssh/mykey2.rsa root@$IPMASTER"
	    echo
	else
	    TAG=$(awk 'FNR == 1 {print $6}' .Info_instances.txt)
    	    echo
    	    echo "  ERROR: Attempt to reconfigure a cluster without a $FILE; check for status; tag is $TAG"
    	    echo
	    exit
	fi
    fi
	
    if [ "$1" = "terminate" ];then
	## Terminate cluster
	FILE=.Info_instances.txt
	if test -f "$FILE"; then
	    echo
	    echo "  Terminating all instances stored in $FILE"
	    cat .Info_instances.txt | xargs -l -n7 bash -c 'aws ec2 terminate-instances --instance-ids $1'  > /dev/null
	    echo "  Removing known hosts stored in $FILE"
	    cat .Info_instances.txt | xargs -l -n7 bash -c 'ssh-keygen -f "/home/tgst/.ssh/known_hosts" -R "$3"'
	    rm $FILE
	    rm .hosts_info.txt
	    rm .hosts_file.txt
	    rm .ip_master.txt
	    echo
	else
	    echo
	    echo "  ERROR: no $FILE with information; need to delete instances manually"
	    echo	
	fi
    fi
    
else
    ## ERROR message
    echo
    echo "  ERROR: Input valid arguments"
    echo "  Usage: Cluster_script (start-spot | start-ondemand | reconfigure | terminate) [omp | acc | mpi]"
    echo
    exit
fi

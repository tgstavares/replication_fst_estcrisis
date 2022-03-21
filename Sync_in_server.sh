#!/bin/bash

SCRIPT_COMMAND='run.sh'
ROOT_DIRECTORY="/home/tgst/Desktop/Estcrisis01"
KEY='/home/tgst/.ssh/mykey2.rsa'

File=.ip_master.txt

if test -f "$File"; then
    IPMASTER=$(cat .ip_master.txt)
else
    echo
    echo "  ERROR: no ip file with information."
    echo	
    exit
fi

if [ "$1" == "mulher" ]; then
    echo
    perl -pe "system 'sleep .04'" ~/.ascii/art/mulher8.txt | lolcat
    echo
fi

if [ "$1" == "run" ]; then
    echo
    echo "  ####> Transfering files to cluster..."
    echo
    cd $ROOT_DIRECTORY
    
    rsync -avz -e "ssh -o StrictHostKeyChecking=no -i $KEY" $ROOT_DIRECTORY/* root@$IPMASTER:/home/ubuntu/works
    echo
    echo "  ####> Compile sources and run binaries..." | lolcat
    echo
    
    ssh -o StrictHostKeyChecking=no -i $KEY root@$IPMASTER "tail -1 ~/.profile > temp && source temp && rm temp && cd /home/ubuntu/works && ./$SCRIPT_COMMAND" | lolcat
    
    echo "  ####> Transfering files from cluster..."
    echo
    rsync -avz -e "ssh -o StrictHostKeyChecking=no -i $KEY" root@$IPMASTER:/home/ubuntu/works/* $ROOT_DIRECTORY/
    echo
fi

if [ "$1" == "cleanup" ]; then
    echo
    echo "  ####> Cleaning up all in cluster..."
    echo
    ssh -o StrictHostKeyChecking=no -i $KEY root@$IPMASTER "rm -r /home/ubuntu/works" | lolcat
fi

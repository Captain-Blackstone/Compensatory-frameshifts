This repository contains scripts I considered to be relevant to understand the algorithm I used in the article and possibly reproduce the results.

Because reproducibility is very important!

However, I highly doubt anyone will go to such trouble just to understand the algorithm... but whatever. If you read this text, you're a good candidate) 

Feel free to contact me through my e-mail.

bibadima@rambler.ru or dmitriy.biba@gmail.com (I really hope to switch to the second adress one day, but currently use the first one).

The script was not intended to be a command line tool. Rather it was intended to operate only once and to deal with all files in the folder
"/mnt/lustre/Blackstone/pipeline_intermediate_data/for_indel_finder/" (or "/mnt/lustre/Blackstone/insect_pipeline_intermediate_data/for_indel_finder/" in case of insects)
and store the results in the folder
"output_indel_finder/"

so, if you don't have such folders (which you most probably don't) on your computer, change the variables in_dir and out_dir in the code.

Also the algorithm uses a phylogenetic tree,which I store in the variable, you know, TREE. You also need toredefine it using your own tree.

One might also wonder, why these scripts and files they produce contain the word "slow". That's because I tried to design a fast version of this script using multiprocessing module, but for some incomprehensible reason it refused to work with more than 1000 files, so I ended up using the slow version.

Also, the script used for insects is considerably different from the first one. That is because I optimized it in terms of speed and code comprehensibiity. Why didn't I make changes to the initial script? Well, that's because it served it purpose and I hardly beleive anyone will ever read it. So, if you find yourself being interested in it, you can e-mail me something like "Man, your script is awful, can you make it, I don't know, better? I really need to understand what's going on there.". Then I will, though pretty reluctantly (who wants to revise his old script!?) respond to your request.

So, I beleive, that's it, my friend. Good luck and have fun!

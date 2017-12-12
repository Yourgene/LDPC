# LDPC
Matlab function implementing the LDPC soft decoding algorithm, with the log sum product method

This code was written by studying the paper of Sarah Johnson (http://sigpromu.org/sarah/SJohnsonLDPCintro.pdf or in resources/SJohnsonLDPCintro.pdf),
so all the maths explanations of my algorithm can be found there.

## USAGE
To use the function, the only paramater is the received sequence (I call it input_frame). You can put bits, or positive/negative floats (it depends of the channel you want to modelize)

## Default values 
- the decoding matrix is :
H = [1 1 0 1 0 0;
     0 1 1 0 1 0;
     1 0 0 0 1 1;
     0 0 1 1 0 1];
 - the channel is BSC (Binary Symmetric Channel) with a 0.2 crossover probability. In the code, you can change it to a AWGN channel.
 
 These default values allow to test examples 2.5 and 2.6 of the paper of Sarah Johnson.
 
 ## Parametrization
 - H, the decoding matrix, can be modified ; the code will adapt. 
 - The impact of the channel on the a priori LLR can be changed. To do that, refer to the differences between fuctions a_priori_log_likehood and a_priori_log_likehood_AWGN
 
 Other variables can be changed too, like the maximum amout of iterations, etc.
 
 
 ## Examples
 
 1. BS Channel with the default H matrix
 command: ldpcdecodersoft([1 0 1 0 1 1]).
 the 1st bit was flipped, and is wrong
 
 2. AWGN Channel with the default H matrix
 You have to replace the call of the function a_priori_log_likehood by a_priori_log_likehood-AWGN.
 command : ldpcdecodersoft([-0.1 0.5 -0.8 1.0 -0.7 0.5]).
 the recommended SBR is 1.25
 
 These 2 examples are from Sarah's paper (examples 2.5 and 2.6)
 

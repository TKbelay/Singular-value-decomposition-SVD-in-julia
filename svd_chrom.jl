# This does SVD of a defined genomic dataset, the dataset should only include genomic data (no header or id data) 
# Run as: julia SVD_chrom.jl <coresnpfile> <allsnpfile> <chromosome to analyse> <fraction wihtin-chr genomic variance to be explained>
infile1=ARGS[1]; #Name of core marker set for the given chromosome
chrom=ARGS[3]; #Chromosome number
infile2=ARGS[2]; #Name of total marker set for the given chromosome

fr_exp=parse(Float64,ARGS[4]); # Fraction of genomic variance per chromosome required to be explained by the chosen components
m0=readdlm("$infile1.$chrom"); # Read in marker data from chr i from file <infile>.<chrom>
m_all0=readdlm("$infile2.$chrom"); # Read in marker data from chr i for all animals

coresize,loci=size(m0); #Define size of core sample population
N=size(m_all0,1); #no. of animals in the total genomic data set for the chromosome

P=0.5*mean(m0,1); # Row vector of allele frequencies in the core
m=m0.-(2*P); # Compute centered marker matrix
rho=2*sum(P.*(1-P));
U,s,V = svd(m); 
wt=(s.^2).*(1/sum(s.^2)); # Fraction of variance explained by each component
yy=0.0;
comps=0;
for j=1:coresize;
	comps=j; 
	yy=yy + wt[j];
	if yy>fr_exp;
	   break;
	end;
end;
# The variable "comps" is the number of chosen components for the chromosome

m_all=m_all0.-(2*P); # Compute centered marker matrix
rho=2*sum(P.*(1-P)); #The part of the VanRaden 1, denominator coming from this chromosome
C=m_all*V[:,1:comps]; #Compute the reduced-dimension score matrix for chromosome i including core and non-core animals

writedlm("scores.txt.$chrom",C,"\t");
writedlm("rho.txt.$chrom", rho,"\t");

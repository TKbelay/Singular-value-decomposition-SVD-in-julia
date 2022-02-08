# This does SVD of a defined genomic dataset, the dataset should only include genomic data (no header or id) 
# Run as: julia SVD_chrom.jl <coresnpfile> <allsnpfile> <chromosome to analyse> <fraction wihtin-chr genomic variance to be explained> (if want run on linux command line)
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

C1=readdlm("score.txt.1");
C2=readdlm("score.txt.2");
C3=readdlm("score.txt.3");
C4=readdlm("score.txt.4");
C5=readdlm("score.txt.5");
C6=readdlm("score.txt.6");
C7=readdlm("score.txt.7");
C8=readdlm("score.txt.8");
C9=readdlm("score.txt.9");
C10=readdlm("score.txt.10");
C11=readdlm("score.txt.11");
C12=readdlm("score.txt.12");
C13=readdlm("score.txt.13");
C14=readdlm("score.txt.14");
C15=readdlm("score.txt.15");
C16=readdlm("score.txt.16");
C17=readdlm("score.txt.17");
C18=readdlm("score.txt.18");
C19=readdlm("score.txt.19");
C20=readdlm("score.txt.20");
C21=readdlm("score.txt.21");
C22=readdlm("score.txt.22");
C23=readdlm("score.txt.23");
C24=readdlm("score.txt.24");
C25=readdlm("core.txt.25");
C26=readdlm("score.txt.26");
C27=readdlm("score.txt.27");
C28=readdlm("score.txt.28");
C29=readdlm("score.txt.29");

#combine components from all chromosomes into none file
C=[C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 C29];

writedlm("pooled_score.txt",C,"\t");

#SVD on pooled components
U,S,V=svd(C);#perform dimension reduction on pooled components (but not further reduction in rank)
T=U*diagm(S);#define components 
s=diagm(S)*diagm(S);#compute square of the singular value (S) i.e S2=T'T

#predict effect of components (aka SNP effects)
IT=Matrix(1.0*I, n,n);#identity matrix with dimension of pooled components(n)
y=readdlm("file.dat");
#t=(inv([s + (h*IT)])*(transpose(T)*y);#as T'T=S2, h=lambda=residual variance/genetic variance
e=h*IT;
f=s+e;
d=inv(f);
t=d*(transpose(T)*y);#predict component effects(as SNP effect(b)=V*t and t=V'b)
g=T*t;#predict individual genetic effects or breeding values

writedlm("ebv.txt",g,"\t")

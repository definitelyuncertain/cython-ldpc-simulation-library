# cimport the Cython declarations for numpy
import numpy as np
from numpy import random as rn
cimport numpy as np
from collections import defaultdict

np.import_array()

def _groupinds(a) :
	D=defaultdict(list)
	for i,val in enumerate(sorted(a)) :
		D[val].append(i)
	return ([min(D[k]) for k in D.keys()],[max(D[k]) for k in sorted(D.keys())])

def _argsort(a) :
	return [a[0] for a in sorted(enumerate(a), key=lambda x:x[1])]

# cdefine the signature of our c function
cdef extern from "LDPC.h":
	void raw_decode_bec(int* bsi, int* bei, int* bls, int* csi, int* cei,
	int n, int nk, int* edg, int* out, int max_itr)
	void printarray(int *a, int n)
	void squarearray(int *a, int n)

# create the wrapper code, with numpy type annotations
#def sigmoid(np.ndarray[double, ndim=1, mode="c"] in_array not None,
#                     np.ndarray[double, ndim=1, mode="c"] out_array not None):
#    sigmoid_array(<double*> np.PyArray_DATA(in_array),
#                <double*> np.PyArray_DATA(out_array),
#                in_array.shape[0])

#def dummywrap(num) :
#	cdef np.ndarray a=np.zeros((1),dtype=int)
#	raw_decode_bec(<int*> np.PyArray_DATA(a),<int*> np.PyArray_DATA(a),\
#	<int*> np.PyArray_DATA(a),<int*> np.PyArray_DATA(a),<int*> np.PyArray_DATA(a),\
#	num,<int*> np.PyArray_DATA(a),<int*> np.PyArray_DATA(a),<int*> np.PyArray_DATA(a))

class LDPC :
	def __init__(self,source) :
		self.row_ind,self.col_ind=[],[]
		if type(source) is str :
			rc=0
			with open(source, 'r') as f:
				dat=[l for l in f]
				self.n,self.nk=np.intc(dat[0].split()[0]),np.intc(dat[0].split()[1])
				for line in dat[1:] :
					self.col_ind=self.col_ind+[int(x) for x in line.split()]
					self.row_ind=self.row_ind+([rc]*len(line.split()))
					rc=rc+1
				f.close()
			(bmin,bmax)=_groupinds(self.col_ind)
			(cmin,cmax)=_groupinds(self.row_ind)
			self.bls=np.array(_argsort(self.col_ind),dtype=np.intc)
			self.bei=np.array(bmax,dtype=np.intc)
			self.bsi=np.array(bmin,dtype=np.intc)
			self.csi=np.array(cmin,dtype=np.intc)
			self.cei=np.array(cmax,dtype=np.intc)
		elif isinstance(source,LDPC) :
			self.n=source.n
			self.nk=source.nk
			self.row_ind=source.row_ind[:]
			self.col_ind=source.col_ind[:]
			self.bls=source.bls.copy()
			self.bei=source.bei.copy()
			self.cei=source.cei.copy()
			self.bsi=source.bsi.copy()
			self.csi=source.csi.copy()
	# Wrapper for BEC decoder routine
	def decodeBEC(self,np.ndarray[int, ndim=1, mode="c"] rxv not None,int max_itr) :
		cdef np.ndarray[int, ndim=1, mode="c"] out = rxv.copy()
		cdef np.ndarray[int, ndim=1, mode="c"] edg = np.zeros(self.bls.size,dtype=np.intc)
		cdef int * outptr = <int*> np.PyArray_DATA(out)
		cdef int * bsiptr = <int*> np.PyArray_DATA(self.bsi)
		cdef int * beiptr = <int*> np.PyArray_DATA(self.bei)
		cdef int * csiptr = <int*> np.PyArray_DATA(self.csi)
		cdef int * ceiptr = <int*> np.PyArray_DATA(self.cei)
		cdef int * blsptr = <int*> np.PyArray_DATA(self.bls)
		cdef int * edgptr = <int*> np.PyArray_DATA(edg)
		raw_decode_bec(bsiptr,beiptr,blsptr,csiptr,ceiptr,self.n,self.nk,edgptr,outptr,max_itr)
		print 'edg',edg
		return out
	# Zero codeword simulation for BEC
	def simBEC(self, double ep, int max_itr, int ber_limit, int fer_limit, int max_times,int diag_interval) :
		cdef np.ndarray[int, ndim=1, mode="c"] out = np.zeros(self.n,dtype=np.intc)
		cdef np.ndarray[int, ndim=1, mode="c"] tosses = np.zeros(self.n,dtype=np.intc)
		cdef np.ndarray[int, ndim=1, mode="c"] edg = np.zeros(self.bls.size,dtype=np.intc)
		# Define pointers to the various arrays
		cdef int * outptr = <int*> np.PyArray_DATA(out)
		cdef int * bsiptr = <int*> np.PyArray_DATA(self.bsi)
		cdef int * beiptr = <int*> np.PyArray_DATA(self.bei)
		cdef int * csiptr = <int*> np.PyArray_DATA(self.csi)
		cdef int * ceiptr = <int*> np.PyArray_DATA(self.cei)
		cdef int * blsptr = <int*> np.PyArray_DATA(self.bls)
		cdef int * edgptr = <int*> np.PyArray_DATA(edg)
		cdef int i=1
		# Storing simulation results
		cdef int bers=0
		cdef int fers=0
		cdef int beritr=0
		while (i<=max_times and bers<ber_limit and fers<fer_limit) :
			# Run bernoulli trials to determine erasure locations
			tosses=rn.binomial(1,1.0-ep,self.n).astype(np.intc)
			out.fill(1)
			edg.fill(0)
			# Erase locations that got 0 in bernoulli trials
			out*=tosses
			# Run decoder
			raw_decode_bec(bsiptr,beiptr,blsptr,csiptr,ceiptr,self.n,self.nk,edgptr,outptr,max_itr)
			beritr=np.where(out==0)[0].size
			bers+=beritr
			if beritr :
				fers+=1
			i=i+1
			if (i%diag_interval)==0 :
				print '''Iteration %d BERs: %d FERs: %d\
				BERate %f FERate %f'''%(i,bers,fers,100*float(bers)/float(self.n*i)\
				,100*float(fers)/float(i))
		return (i,bers,fers,100*float(bers)/float(self.n*i)\
		,100*float(fers)/float(i))
# End of class definition
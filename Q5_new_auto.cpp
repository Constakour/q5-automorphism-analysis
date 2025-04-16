#include <bits/stdc++.h>
using namespace std;

/*
  Enumerates all pairs-only half-sets in Q_5:
    - Q_5 has 32 vertices, labeled 0..31 (5-bit).
    - There are 16 antipodal pairs => (v, v^31).
    - We'll skip the pair containing vertex31, so we choose 8 from the other 15 => 6435 combos.

  For each half-set S (16 vertices):
    1) If a single XOR flip can do S->~S => skip
    2) Else, test "flip->perm": Tset = perm( flip(S) ).
       If Tset==~S, we store S in q5_required_permutation.json with the found transformation.

  This code ensures only bits 0..31 matter, because we build sets in a 32-bit domain. 
  The final JSON contains:
    {
      "mask_hex": "0x????????",
      "vertices": [... 16 integers ...],
      "flip_mask": <0..31>,
      "perm": [p0,p1,p2,p3,p4]
    }
  so that later we can do the same "flip->perm" in Python or Sage.

  Note: The "flip->perm" means we do X=flipMask XOR on each vertex in S, then apply a coordinate permutation perm[] to reorder bits.

  We'll skip writing sets that only require a single translation (XOR).
*/

static constexpr int N = 32; // 2^5 for Q_5

struct Pair {
    int a, b;
};

// Apply an XOR flip
unsigned int applyXor(unsigned int S, unsigned int flipMask){
    unsigned int T=0U;
    for(int i=0;i<N;i++){
        if(S & (1U<<i)){
            int x = i ^ flipMask; // guaranteed <32 if i<32 & flipMask<32
            T |= (1U << x);
        }
    }
    return T;
}

// Permute bits: perm[i] = oldBit => newBit i
unsigned int applyPermutation(unsigned int S, const array<int,5> &p){
    unsigned int P=0U;
    for(int v=0; v<N; v++){
        if(S & (1U<<v)){
            int newV=0;
            for(int b=0;b<5;b++){
                int oldBit = (v >> p[b]) & 1;
                newV |= (oldBit << b);
            }
            // newV < 32
            P |= (1U << newV);
        }
    }
    return P;
}

// single 5-bit translation check
bool checkSingleTranslation(unsigned int S){
    unsigned int comp = (~S) & 0xFFFFFFFFU; 
    for(unsigned int flip=0; flip<32; flip++){
        unsigned int T= applyXor(S, flip);
        if(T == comp) return true;
    }
    return false;
}

// We'll define permutations => "flip->perm" approach
static vector< array<int,5> > allPerms;
static bool builtPerms=false;

bool checkFlipPerm(unsigned int S,
                   unsigned int &out_flip,
                   array<int,5> &out_perm)
{
    unsigned int comp = (~S) & 0xFFFFFFFFU;
    if(!builtPerms){
        vector<int> base{0,1,2,3,4};
        sort(base.begin(), base.end());
        do {
            array<int,5> p;
            for(int i=0;i<5;i++){
                p[i] = base[i];
            }
            allPerms.push_back(p);
        } while(next_permutation(base.begin(), base.end()));
        builtPerms=true;
    }

    // for each flip<32, for each perm => if perm(flip(S))==comp => success
    for(auto &p : allPerms){
        for(unsigned int flip=0; flip<32; flip++){
            unsigned int X= applyXor(S, flip);
            unsigned int T= applyPermutation(X, p);
            if(T == comp){
                out_flip= flip;
                out_perm= p;
                return true;
            }
        }
    }
    return false;
}

//////////////////////////////////////////////////
// We'll build the 16 antipodal pairs for Q_5
// Then skip the pair containing vertex31 => pick 8 from the other 15 => 6435 sets
//////////////////////////////////////////////////

// next combination with same popcount => Gosper
unsigned long long nextCombPopcount(unsigned long long x){
    if(x==0ULL) return 0ULL;
    unsigned long long u= x & -x;
    unsigned long long v= u + x;
    if(v==0ULL) return 0ULL;
    x= v + (((v ^ x)/u)>>2);
    return x;
}

// binomial(15,8)=6435
unsigned long long binom(int n,int k){
    long double r=1.0;
    for(int i=1; i<=k; i++){
        r*= (n-k+i)/(long double)i;
    }
    return (unsigned long long) llroundl(r);
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // build pairs
    vector<bool> used(N,false);
    vector<Pair> pairsList;
    for(int i=0;i<N;i++){
        if(!used[i]){
            int j= i ^ 31; // 31=11111b => antipode in Q_5
            used[i]=true; used[j]=true;
            pairsList.push_back({ min(i,j), max(i,j) });
        }
    }
    if(pairsList.size()!=16){
        cerr<<"Error building pairs??\n";
        return 1;
    }
    // find pair with vertex31 => skip it
    int fixPairIdx=-1;
    for(int i=0;i<16;i++){
        if(pairsList[i].b==31){
            fixPairIdx=i;
            break;
        }
    }
    if(fixPairIdx<0){
        cerr<<"No pair containing vertex31??\n";
        return 1;
    }
    // build a mask for each pair
    vector<unsigned int> pairMask(16,0U);
    for(int i=0;i<16;i++){
        unsigned int mm= (1U<< pairsList[i].a) | (1U<< pairsList[i].b);
        pairMask[i]= mm;
    }

    // we pick 8 from the other 15 => 6435 combos
    unsigned long long total= binom(15,8);
    cout<<"Enumerating "<< total <<" distinct half-sets in Q_5 (mod complement).\n";

    // define idxList => skipping fixPairIdx
    vector<int> idxList; idxList.reserve(15);
    for(int i=0;i<16;i++){
        if(i!= fixPairIdx) idxList.push_back(i);
    }

    auto buildSetMask = [&](unsigned long long c){
        unsigned int S=0U;
        while(c){
            unsigned long long v= (c & -c);
            int b= __builtin_ctzll(c);
            c^= v;
            int pairIdx= idxList[b];
            S |= pairMask[pairIdx];
        }
        return S;
    };

    unsigned long long startC= (1ULL<<8)-1ULL; //lowest => 0xFF
    unsigned long long endC= startC<<(15-8);

    ofstream fout("q5_required_permutation.json");
    if(!fout){
        cerr<<"Cannot open q5_required_permutation.json\n";
        return 1;
    }
    fout<<"{\n  \"sets\": [\n";

    bool firstOutput=true;
    unsigned long long enumerated=0ULL;
    unsigned long long comb= startC;
    const unsigned long long PRINT_INTERVAL=500ULL;

    while(true){
        enumerated++;
        unsigned int S= buildSetMask(comb); // a 16-vertex half-set in Q_5

        // if single translation => skip
        if(!checkSingleTranslation(S)){
            // else do flip->perm
            unsigned int found_flip=0;
            array<int,5> found_perm;
            if(checkFlipPerm(S, found_flip, found_perm)){
                // store in JSON
                if(!firstOutput) fout<<",\n";
                firstOutput=false;

                // gather vertices
                vector<int> verts;
                for(int v=0; v<N; v++){
                    if(S & (1U<<v)) verts.push_back(v);
                }
                // build mask in hex
                char buf[16];
                sprintf(buf, "0x%08X", S);

                fout<<"    {\n";
                fout<<"      \"mask_hex\": \""<< buf <<"\",\n";
                fout<<"      \"vertices\": [";
                for(size_t i=0;i<verts.size();i++){
                    if(i>0) fout<<", ";
                    fout<< verts[i];
                }
                fout<<"],\n";
                fout<<"      \"flip_mask\": "<< found_flip <<",\n";
                fout<<"      \"perm\": [";
                for(int i=0;i<5;i++){
                    if(i>0) fout<<", ";
                    fout<< found_perm[i];
                }
                fout<<"]\n    }";
            }
        }

        if((enumerated%PRINT_INTERVAL)==0ULL){
            cout<<"Processed "<<enumerated<<"/"<<total<<" sets...\n";
        }
        if(comb== endC) break;
        unsigned long long nxt= nextCombPopcount(comb);
        if(nxt==0ULL) break;
        comb= nxt;
    }

    fout<<"\n  ]\n}\n";
    fout.close();

    cout<<"Done. Processed "<< enumerated <<" half-sets.\n";
    cout<<"Results stored in q5_required_permutation.json\n";
    return 0;
}

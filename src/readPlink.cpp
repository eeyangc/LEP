
#include "readPlink.hpp"



float normalCFD(float value)
{
    return 0.5 * erfc(-value * M_SQRT1_2);
}


/*get the positions of the identifiers in the column names*/
Col<int> getPositions(vector <string> fields, vector<string>& identifiers){

    uword size = identifiers.size();
    Col<int> pos(size);
    for(uword i = 0; i < size; i++){
        pos[i] = -1;
        for(uword j = 0; j < fields.size(); j++){
            if(fields[j].compare(identifiers[i]) == 0){
                pos[i] = (int)j ;
                break;
            }
        }

    }
    return pos;

}

GenoInfo::GenoInfo(string stringname) {

    string famfile = stringname;
    famfile += ".fam";
    string bimfile = stringname;
    bimfile += ".bim";

    this -> N =  getLineNum(famfile);
    this -> P =  getLineNum(bimfile);

    clock_t t1 = clock();
    Chroms chroms(bimfile, P);
    int phenotype_pos = 5;
    vec y = read_phenotypes(famfile, N, phenotype_pos);
    long long size = (long long)N * (long long)P;
    unsigned* X = new unsigned[ size];
    t1 = clock();
    readPlink(stringname,N, P, X);
    arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
    this -> X = *Xdata;
    this -> X.replace(3, 0);
    this -> chroms = chroms;
    this -> y = y;
    cout<<"Finish Reading Plink file, time elapsed:"<< (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
    cout <<"Sample Size = " <<  N << " SNP Number = " << P << endl;
    delete[] X;

}




int snps_overlap(vector<SNP>& chrom_x_i, vector<SNP>& chrom_s_i, Col<uword>& xindex, Col<uword>& sindex){

    vector<SNP> common_snp_in_x;
    vector<SNP> common_snp_in_ss;
    sort(chrom_s_i.begin(), chrom_s_i.end());
    sort(chrom_x_i.begin(), chrom_x_i.end());

    //  cout << chrom_x_i[0] << "  " << chrom_x_i[1] << endl;
    set_intersection(chrom_x_i.begin(),chrom_x_i.end(),chrom_s_i.begin(),chrom_s_i.end(),back_inserter(common_snp_in_x));

    set_intersection(chrom_s_i.begin(),chrom_s_i.end(), chrom_x_i.begin(),chrom_x_i.end(),back_inserter(common_snp_in_ss));

    //  xindex.resize(common_snp_in_x.size());
    Mat<uword> indexPair(common_snp_in_x.size(), 2);
    for(int i = 0; i < common_snp_in_x.size(); i++){
        indexPair(i,0) = (uword)common_snp_in_x[i].idx;
        indexPair(i,1) = (uword)common_snp_in_ss[i].idx;
    }

    uvec sorted_index = sort_index(indexPair.col(0));
    indexPair = indexPair.rows(sorted_index);
    xindex = indexPair.col(0);
    sindex = indexPair.col(1);

    return (int)common_snp_in_x.size();

}

vec read_phenotypes(string filename, int N, int phenopos){
    std::ifstream ifs(filename.c_str());

    std::string line;
    vec y(N);
    string snpname;
    float pheno = 0;
    string tmp_str;
    int idx = 0;
    while(std::getline(ifs, line)) // read one line from ifs
    {
        std::istringstream iss(line); // access line as a stream
        vector<string> fields;
        boost::split( fields, line, boost::is_any_of(" \t"));
        pheno = (double)atof(fields[phenopos].c_str());
        y[idx] = pheno;
        idx++;
    }
    ifs.close();
    return y;
}

Chroms read_snpnames(string filename, int P){
  std::ifstream ifs(filename.c_str());

  std::string line;
  Chroms chroms(P);//snps(P);
  int chromsome;
  string snpname;
  vector <string> fields;
  int i = 0;
  while(std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    boost::split( fields, line, boost::is_any_of(" \t *"));
    chromsome = (int)atoi(fields[0].c_str());
    snpname = fields[1];
    int a1 = *fields[4].c_str();
    int a2 = *fields[5].c_str();
    chroms.chromsome[i] = chromsome;
    SNP snp(snpname, (int)i, from_x);
    chroms.snps.push_back(snp);
    chroms.A1Array[i] = a1;
    chroms.A2Array[i] = a2;
    i++;

  }
  ifs.close();
  return chroms;
}

int chroms_overlap(Chroms& chrom_x, Chroms& chrom_s,Col<uword>& xindex, Col<uword>& sindex){
  int total_overlap = 0;
  total_overlap = snps_overlap(chrom_x.snps, chrom_s.snps,  xindex, sindex);
  return total_overlap;

}

#define keep_cols 0
#define keep_rows 1

template <typename T>
void keep_indices(arma::Mat<T>& mat_obj, int type, uvec index)  {
    arma::Mat<T> mat_tmp;
    if(type == keep_cols){
       mat_tmp = mat_obj.cols(index);
    }else{
       mat_tmp = mat_obj.rows(index);
    }
    mat_obj.reset();
    mat_obj = mat_tmp;
}

void Summary::cal_overlap(GenoInfo& genoinfo)
{
    Col<uword> xindex;
    Col<uword> sindex;

    Col<int> direction;
    direction.reset();

    chroms_overlap(genoinfo.chroms, *this -> chroms, xindex, sindex);
    keep_indices(genoinfo.X, keep_cols, xindex);
    keep_indices((*this -> lpsummary), keep_rows, sindex);


    Chroms  chrom((int)xindex.size());
    chrom.chromsome = genoinfo.chroms.chromsome.elem(xindex);
    chrom.A1Array = genoinfo.chroms.A1Array.elem(xindex);
    chrom.A2Array = genoinfo.chroms.A2Array.elem(xindex);
    for(int i = 0; i < chrom.P; i++){
        chrom.snps[i] = genoinfo.chroms.snps[xindex[i]];
    }
    genoinfo.chroms.clear();
    genoinfo.chroms = chrom;
    genoinfo.xindex = xindex;
    genoinfo.N = (int)genoinfo.X.n_rows;
    genoinfo.P = (int)genoinfo.X.n_cols;
}


Summary::Summary(vector<string>& snps, mat* lp_summary){
  int type = beta_ss;
  this -> P = snps.size();
  this -> chroms = new Chroms(P);
  this -> type = type;
  this ->lpsummary = lp_summary;
  for(uword k = 0; k < P; k++){
    SNP snp(snps[k], k, from_ss);
    this -> chroms -> snps.push_back(snp);
  }
  this -> chrom_no = 1;
}

bool SNP::operator<(const SNP& obj) const{
    return this -> name < obj.name;
}
bool SNP::operator>(const SNP& obj) const{
    return this -> name > obj.name;
}

bool SNP::operator != (const SNP& obj) const{
    return this -> name.compare(obj.name) > 0;
}


bool SNP::operator == (const SNP& obj) const{
    return this -> name.compare(obj.name) == 0;
}




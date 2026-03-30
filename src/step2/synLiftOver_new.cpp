#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <regex>
#include <dirent.h>
#include <unistd.h>
#include <list>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <limits.h>

using namespace std;


vector < string > SPLIT(string s, char delimiter = '1') {
    vector < string > return_arr;
    istringstream iss(s);
    if (delimiter == '1') { //default, split with spaces
        do {
            string sub;
            iss >> sub;
            return_arr.push_back(sub);
        } while (iss);

        return_arr.pop_back();
        return return_arr;
    } else {
        string temp;
        while (getline(iss, temp, delimiter)) {
            return_arr.push_back(temp);
        }
        return return_arr;
    }
}
int STR2INT(string s) {
    return atoi(s.c_str());
}
string SLICING (vector <string> input, int restStartingIdx){
    string rest_str;
    int size = input.size();
    for(int i = restStartingIdx ; i < size ; i ++){
        rest_str += input[i]+"\t";
    }

    return rest_str;
}

// 
struct query {
    string chr;
    int start;
    int end;
    string rest_information;
    struct axt_list *stAxt = nullptr;
    struct axt_list *edAxt = nullptr;
};
struct axt_list {
    string ref_seq;
    string tar_seq;
    int refStart;
    int refEnd;
    string tarChr;
    int tarStart;
    int tarEnd;
    int thisIdx;
    string tarDir;
    struct axt_list *next =nullptr;
};

// 
vector < query > readQuery(string& path){
    struct stat sb{};
    string res;
    vector < query > return_list;
    vector <query> thischr_tmp;

    FILE* input_file = fopen(path.c_str(), "r");
    if (input_file == nullptr) {
        perror("fopen");
    }

    stat(path.c_str(), &sb);
    res.resize(sb.st_size);
    fread(const_cast<char*>(res.data()), sb.st_size, 1, input_file);
    fclose(input_file);

    vector <string> all_query = SPLIT(res,'\n');
    vector < string > tmp_arr;
    string lastchr = "NULL";
    int all_query_lines = all_query.size();
    for (int i = 0 ; i < all_query_lines ; i+=1){ 
        tmp_arr = SPLIT(all_query[i]);
        struct query cur_query;
        cur_query.chr = tmp_arr[0];
        cur_query.start = STR2INT(tmp_arr[1]);
        cur_query.end = STR2INT(tmp_arr[2]);
        cur_query.rest_information = SLICING(tmp_arr,3);
        return_list.push_back(cur_query);

    }  
    return return_list;
}

map < string, int > readSize(string& path){
    struct stat sb{};
    string res;
    map < string, int > return_map;

    FILE* input_file = fopen(path.c_str(), "r");
    if (input_file == nullptr) {
        perror("fopen");
    }

    stat(path.c_str(), &sb);
    res.resize(sb.st_size);
    fread(const_cast<char*>(res.data()), sb.st_size, 1, input_file);
    fclose(input_file);

    vector <string> all_Size = SPLIT(res,'\n');
    vector <string > tmp_arr;
    int all_size_lines = all_Size.size();
    for (int i = 0 ; i < all_size_lines ; i+=1){
        tmp_arr=SPLIT(all_Size[i]);
        return_map[tmp_arr[0]] = STR2INT(tmp_arr[1]);
    }
    return return_map;
}
struct axt_list * readAxt(string path) { // axt should not have any blank lines 
    struct axt_list * return_axt = new axt_list();
    struct axt_list * last = new axt_list();

    int first_flag = 1;
    string line;
    ifstream input_file(path);
    int countF = 0;

    vector < string > tmp_arr;
    struct axt_list *curAxt = nullptr;
    while (getline(input_file, line)){
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (countF ==0){
            struct axt_list *thisAxt = new axt_list();
            curAxt = thisAxt;
            tmp_arr = SPLIT(line);
            thisAxt->thisIdx = STR2INT(tmp_arr[0]);
            thisAxt->refStart = STR2INT(tmp_arr[2]);
            thisAxt->refEnd = STR2INT(tmp_arr[3]);
            thisAxt->tarChr = tmp_arr[4];
            thisAxt->tarStart = STR2INT(tmp_arr[5]);
            thisAxt->tarEnd = STR2INT(tmp_arr[6]);
            thisAxt->tarDir = tmp_arr[7];
            thisAxt->next = NULL;
        }
        if (countF==1){
            curAxt->ref_seq = line;
        }
        if (countF == 2){
            curAxt->tar_seq = line;
            if (first_flag==1){
                return_axt = last = curAxt;
                first_flag = 0;
            }else{
                last->next = curAxt;
                last = curAxt;
            }
            countF = 0;
        }else{
            countF ++;
        }

        
    }
    input_file.close();
    return return_axt;
}
void convertCdnt (string *refseq, string *tarseq, int refst, int tarst, int *liftTar1, int *liftTar2){
    refst -= 1;
    tarst -= 1;
    int done1 = (*liftTar1);
    int done2 = (*liftTar2);
    int converted1 = 0;
    int converted2 = 0;

    if ((*refseq).length() != (*tarseq).length()){
        cout <<"err" <<endl;
    }
    string gap = "-";
    for(int idx = 0 ; idx <= (*refseq).length() ; idx ++ ){
        if (done1 == refst){
            converted1 = tarst;
            done1 = 1;
        }
        if (done2 == refst){
            converted2 = tarst;
            done2 = 1;
        }
        if (done1 == 1 && done2 ==1){
            (*liftTar1) = converted1;
            (*liftTar2) = converted2;
            break;
        }

        if ((*refseq)[idx] == '-'){
            tarst ++;
        }else if((*tarseq)[idx] == '-'){
            refst ++;
        }else{
            refst ++;
            tarst ++;
        }
    }
    return;
}

void convertDir(int *liftTar1, int *liftTar2, int chr_size){
    int save_1 = (*liftTar1);
    int save_2 = (*liftTar2);
    (*liftTar1) = chr_size - save_2;
    (*liftTar2) = chr_size - save_1;
    return;
}
// 
int caseClassifier(int ref_s, int ref_e, int q_s, int q_e, float threshold){
    ref_s -= 1;
    int Ccase = 0;
    int qlen = float(q_e - q_s);
    int ovlpLen = 0;
    float ovlpratio = 1;

    if (ref_s <= q_s && q_e <= ref_e){
        Ccase=1;
    }else if(ref_s >= q_s && ref_e <= q_e){
        ovlpLen = ref_e - ref_s;
        ovlpratio = float(ovlpLen) / qlen;
        if(ovlpratio > threshold){
            Ccase = 4;
        }
    }else if(ref_s <= q_s && ref_e >= q_s){
        if(ref_s <= q_e && ref_e <= q_e){
            ovlpLen = ref_e - q_s;
            ovlpratio = float(ovlpLen)/qlen;
            if(ovlpratio > threshold){
                Ccase = 2;
            }
        }
    }else if(ref_s >= q_s && ref_e >= q_s){
        if(ref_s <= q_e && ref_e >= q_e){
            ovlpLen = q_e - ref_s;
            ovlpratio = float(ovlpLen) / qlen;
            if(ovlpratio > threshold){
                Ccase=3;
            }
        }
    }
    return Ccase;
}
int processCase2(const vector<axt_list*> &stCandidates,const vector<axt_list*> &edCandidates,const query &q,map<string,int> &chr_size, int duplication_allow, axt_list* allAxt){
    int q_start = q.start;
    int q_end   = q.end;
    int idx = 0;
    int return_value = 0;
    for (auto stAxt : stCandidates) {
        for (auto edAxt : edCandidates) {
            if (stAxt->thisIdx >= edAxt->thisIdx) {
                continue;
            }
            if (stAxt->tarChr != edAxt->tarChr) {
                continue;
            }
            if (stAxt->tarDir != edAxt->tarDir) {
                continue;
            }
            if (stAxt->refStart > edAxt->refStart) {
                continue;
            }
            if (stAxt->tarStart > edAxt->tarStart) {
                continue;
            }
            bool is_valid_pair = true;
            for (axt_list* other = allAxt; other != nullptr; other = other->next) {
                if (other == stAxt || other == edAxt) continue;
                bool ref_in_between =
                    (other->refStart >= stAxt->refEnd) ||
                    (other->refEnd <= edAxt->refStart);

                bool query_in_between =
                    (other->tarStart >= stAxt->tarEnd) ||
                    (other->tarEnd <= edAxt->tarStart);

                if (ref_in_between || query_in_between) {
                    is_valid_pair = false;
                    break;
                }
            }
            if (!is_valid_pair) {
                continue;
            }

            int lift_st = q_start;
            int lift_ed = q_end; 
            int dummy = 1; 

            convertCdnt(&stAxt->ref_seq, &stAxt->tar_seq,stAxt->refStart, stAxt->tarStart,&lift_st, &dummy);
            dummy = 1; 
            convertCdnt(&edAxt->ref_seq, &edAxt->tar_seq,edAxt->refStart, edAxt->tarStart,&lift_ed, &dummy);
            if (edAxt->tarDir[0] == '-') {
                convertDir(&lift_st, &lift_ed, chr_size[edAxt->tarChr]);
            }
            if (lift_st >= lift_ed) {
                //cerr << "[Case2 warn] converted coords reversed or invalid\n";
                continue;
            }
            cout << edAxt->tarChr << "\t"<< lift_st << "\t"<< lift_ed << "\t"<< q.rest_information <<endl;
            cerr <<q.chr<<"\t"<<q.start<<"\t"<<q.end<<"\t" << edAxt->tarChr<<"\t" <<lift_st << "\t" <<lift_ed <<"\t"<<q.rest_information<<endl;
            return_value = 1;
            if(duplication_allow != 0){
                idx ++;
            }else{
                break;
            }
        }
    }
    return return_value;
}

int processCase2_relaxed(const vector<axt_list*> &stCandidates,const vector<axt_list*> &edCandidates,const query &q,map<string,int> &chr_size, int duplication_allow, axt_list* allAxt){
    int q_start = q.start;
    int q_end   = q.end;
    int idx = 0;
    int return_value = 0;
    for (auto stAxt : stCandidates) {
        for (auto edAxt : edCandidates) {
            if (stAxt->thisIdx >= edAxt->thisIdx) {
                continue;
            }
            if (stAxt->tarChr != edAxt->tarChr) {
                continue;
            }
            if (stAxt->tarDir != edAxt->tarDir) {
                continue;
            }
            if (stAxt->refStart > edAxt->refStart) {
                continue;
            }
            if (stAxt->tarStart > edAxt->tarStart) {
                continue;
            }

            int lift_st = q_start;
            int lift_ed = q_end; 
            int dummy = 1; 

            convertCdnt(&stAxt->ref_seq, &stAxt->tar_seq,stAxt->refStart, stAxt->tarStart,&lift_st, &dummy);
            dummy = 1; 
            convertCdnt(&edAxt->ref_seq, &edAxt->tar_seq,edAxt->refStart, edAxt->tarStart,&lift_ed, &dummy);
            if (edAxt->tarDir[0] == '-') {
                convertDir(&lift_st, &lift_ed, chr_size[edAxt->tarChr]);
            }
            if (lift_st >= lift_ed) {
                //cerr << "[Case2 warn] converted coords reversed or invalid\n";
                continue;
            }
            cout << edAxt->tarChr << "\t"<< lift_st << "\t"<< lift_ed << "\t"<< q.rest_information <<endl;
            cerr <<q.chr<<"\t"<<q.start<<"\t"<<q.end<<"\t" << edAxt->tarChr<<"\t" <<lift_st << "\t" <<lift_ed <<"\t"<<q.rest_information<<endl;
            return_value = 1;
            if(duplication_allow != 0){
                idx ++;
            }else{
                break;
            }
        }
    }
    return return_value;
}




// 
int main(int argc, char * argv[]){
    string cur_exec_name = argv[0];
    vector < string > argList(argv, argv + argc);
    string f_query = argList[1];
    string thischr = argList[2];
    string f_targetSize = argList[3];
    float param_minMatch = std::stof(argList[4]);
    string f_axt = argList[5];
    int duplication_allow = STR2INT(argList[6]);
    int relax = STR2INT(argList[7]);
    
    struct axt_list *thisChrAxt = NULL;
    
    vector < query > input_list;
    vector < query > rest_input_list;
    map < string, int > chr_size;

    input_list = readQuery(f_query);
    chr_size = readSize(f_targetSize);
    thisChrAxt = readAxt(f_axt);
    
    for(int i = 0 ; i < input_list.size(); i++){
        query &q = input_list[i];
        vector<axt_list*> case1_blocks; 
        vector<int> case1_list;

        vector<axt_list*> stCandidates;
        vector<axt_list*> edCandidates;
        for (axt_list* cur = thisChrAxt; cur != NULL; cur = cur->next) {
            int ccase = caseClassifier(cur->refStart, cur->refEnd, q.start, q.end, param_minMatch);
            if(ccase != 0){
                case1_blocks.push_back(cur); 
                case1_list.push_back(ccase); 
            }else{
                if ((cur->refStart -1) <= q.start && q.start <= cur->refEnd) { //0base - 1base comparing
                    stCandidates.push_back(cur);
                }
                if ((cur->refStart-1) <= q.end && q.end <= cur->refEnd) {
                    edCandidates.push_back(cur);
                }
            }
        }
        int case2_processed = 0;
        if(relax == 0){
            case2_processed = processCase2(stCandidates, edCandidates, q, chr_size, duplication_allow, thisChrAxt);
        }else{
            case2_processed = processCase2_relaxed(stCandidates, edCandidates, q, chr_size, duplication_allow, thisChrAxt);
        }
        if (!case1_blocks.empty()) {
            int idx = 0;
            for (auto &block: case1_blocks) {
                int lift_st = q.start;
                int lift_ed = q.end;
                int ccase = case1_list[idx];
                switch (ccase) {
                    case 1:
                        convertCdnt(&block->ref_seq, &block->tar_seq,block->refStart, block->tarStart,&lift_st, &lift_ed);
                        break;
                    case 2:
                        lift_ed = 1;
                        convertCdnt(&block->ref_seq, &block->tar_seq,block->refStart, block->tarStart,&lift_st, &lift_ed);
                        lift_ed = block->tarEnd;
                        break;
                    case 3:
                        lift_st = 1 ;
                        convertCdnt(&block->ref_seq, &block->tar_seq,block->refStart, block->tarStart,&lift_st, &lift_ed);
                        lift_st = block->tarStart -1 ;
                        break;

                    case 4:
                        lift_st = block->tarStart -1 ;
                        lift_ed = block->tarEnd;
                        break;
                }
                if (block->tarDir[0] == '-') {
                    convertDir(&lift_st, &lift_ed, chr_size[block->tarChr]);
                }
                if(!(case2_processed == 1 && duplication_allow == 0)){
                    cout << block->tarChr<<"\t" <<lift_st << "\t" <<lift_ed <<"\t"<<q.rest_information<<"\t"<<endl;
                    cerr <<q.chr<<"\t"<<q.start<<"\t"<<q.end<<"\t" << block->tarChr<<"\t" <<lift_st << "\t" <<lift_ed <<"\t"<<q.rest_information<<endl;
                }
                if(duplication_allow != 0){
                    idx ++;
                }else{
                    break;
                }
            }
        }
        
    }
    return 0;
}

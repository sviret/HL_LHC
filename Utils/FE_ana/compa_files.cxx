// Base class for output comparison
// 
// Output is stored in binary text files corresponding to the front-end data format 

// Meaning of the parameters
//
//
// verilog  : Output of the verilog model (text file)
// simu     : Output of the simu model (text file define at step 7 of
//            http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto620
// 
// outfile  : name of the output text file containing comparison output
// type     : 0 for trigger data, 1 for L1

#include "compa_files.h"

// Main constructor

compa_files::compa_files(std::string verilog, std::string simu,
			 std::string outfile, int type, int MPA)
{
  m_out   = outfile; 
  m_insim = simu; 
  m_indat = verilog; 
  m_isMPA = MPA;

  std::cout << "...Comparing output files " << verilog 
	    << " and " << simu << " containing ";

  (type == 0)
    ? std::cout << "trigger data "
    : std::cout << "L1 raw data ";

  std::cout << std::endl;
      
  (type == 0)
    ? compa_files::compa_TRG()
    : compa_files::compa_L1();
}


//////////////////////////////////////////////////////////////////////
//
// This method reads and decode the concentrator output
// then prepare an output file to compare it to the input data 
// file initially sent
//
// First decode the trigger block, then decode the L1 raw block, either sparsified or not 
//
//
//////////////////////////////////////////////////////////////////////

void compa_files::compa_TRG()
{
  std::cout << "Into compa TRG" << std::endl;

    fver.open(m_indat); // open a file
    fsim.open(m_insim);    // open a file
    
    std::vector<int> sec_ver;
    std::vector<int> sec_sim;
    std::vector<int> lines;

    if (!fver.good() || !fsim.good())
    {
      std::cout << "Please provide valid filenames " << std::endl;
      return; // exit if file not found
    }    

    int line_ver = 0;
    int line_sim = 0;

    char bufver[512];
    char bufsim[512];

    char* tokenver[50]; // initialize to 0
    char* tokensim[50]; // initialize to 0

    int n,m;
    
    int tic = 0;
    int toc = 0;
    
    sec_sim.clear();
    sec_ver.clear();
    
    // read each line of the file
    while (!fver.eof() || !fsim.eof())
    {
        // read an entire line into memory

        if (tic<=toc) fver.getline(bufver, 512);
        if (toc<=tic) fsim.getline(bufsim, 512);

        ++line_ver;
        
        if (line_ver<8) continue;
        
        if (strlen(bufver)<5)
        {
            if (tic<=toc) ++tic;
            if (tic==toc) continue;
        }
            
        if (strlen(bufsim)<5)
        {
            if (toc<=tic) ++toc;
            if (tic==toc) continue;
        }
        
       // cout << line_ver << " / " << strlen(bufver) << " / " << strlen(bufsim) << endl;
        //cout << tic << " / " << toc << endl;
        
        if (tic!=toc) continue;

        
        // array to store memory addresses of the tokens in buf

        // parse the line
        tokenver[0] = strtok(bufver, " |"); // first token
        
        if (tokenver[0] != NULL) // zero if line is blank
        {
            for (n = 1; n < 50; n++)
            {
                tokenver[n] = strtok(NULL, " |"); // subsequent tokens
                if (!tokenver[n]) break; // no more tokens
            }
        }
   
        tokensim[0] = strtok(bufsim, " |"); // first token
        
        if (tokensim[0] != NULL) // zero if line is blank
        {
            for (m = 1; m < 50; m++)
            {
                tokensim[m] = strtok(NULL, " |"); // subsequent tokens
                if (!tokensim[m]) break; // no more tokens
            }
        }
        
        //if (n!=m) break;
        
	//	std::cout << n << "/" << m << std::endl;
        
        if (n==3 && m==5)
        {
            std::string tver(tokenver[1]);
            std::string tsim(tokensim[4]);
            
	    std::cout << tver << " / " << tsim << std::endl;
            for (unsigned int i = 0; i < tver.size(); i++)
            {
                sec_ver.push_back(tver.at(i)-48);
                sec_sim.push_back(tsim.at(i)-48);
                lines.push_back(line_ver);
            }
        }
             
        strcpy(bufver,"");
        strcpy(bufsim,"");
    }
    
    for (unsigned int i = 0; i < sec_ver.size(); i++)
    {
        if (sec_sim[i]!=sec_ver[i])
        {
	  std::cout << lines[i] << "/" << sec_ver[i] << "/" << sec_sim[i] << std::endl;
        }
    }
    
    std::cout << std::endl;

 
    return;
}



void compa_files::compa_L1()
{
  std::cout << "Into compa L1" << std::endl;

  fver.open(m_indat);    // open a file containing verilog info
  fsim.open(m_insim);    // open a file containing simu info
    
  std::vector<int> sec_ver; // Contains the full data flow for one chip
  std::vector<int> sec_sim; //
  std::vector<int> lines;

  if (!fver.good() || !fsim.good())
  {
    std::cout << "Please provide valid filenames " << std::endl;
    return; // exit if file not found
  }    

  int line_ver = 0;
  int line_sim = 0;

  char bufver[512];
  char bufsim[512];

  char* tokenver[50]; // initialize to 0
  char* tokensim[50]; // initialize to 0

  int n,m;
    
  int tic = 0;
  int toc = 0;
    
  sec_sim.clear();
  sec_ver.clear();
    
  int chipn = 0;

  // read each line of the file
  while (!fver.eof() || !fsim.eof())
  {
    // read an entire line into memory
    
    if (tic<=toc) fver.getline(bufver, 512);
    if (toc<=tic) fsim.getline(bufsim, 512);

    ++line_ver;
        
    if (line_ver<8) continue;
        
    if (strlen(bufver)<5) // This thing is to detect when we deal with a new chip
    {
      if (tic<=toc) ++tic;
      if (tic==toc) continue;
    }
            
    if (strlen(bufsim)<5)
    {
      if (toc<=tic) ++toc;
      if (tic==toc) continue;
    }
                
    if (tic!=toc) continue;
        
    if (chipn!=toc) // End of the recording for one chip
    {
      compa_files::do_ana(sec_ver,sec_sim,chipn);
      sec_ver.clear();
      sec_sim.clear();
      ++chipn;
    }
    
    // parse the line
    tokenver[0] = strtok(bufver, " |"); // first token
        
    if (tokenver[0] != NULL) // zero if line is blank
    {
      for (n = 1; n < 50; n++)
      {
	tokenver[n] = strtok(NULL, " |"); // subsequent tokens
	if (!tokenver[n]) break; // no more tokens
      }
    }
   
    tokensim[0] = strtok(bufsim, " |"); // first token
        
    if (tokensim[0] != NULL) // zero if line is blank
    {
      for (m = 1; m < 50; m++)
      {
	tokensim[m] = strtok(NULL, " |"); // subsequent tokens
	if (!tokensim[m]) break; // no more tokens
      }
    }
     
    if (n!=m) continue;
    
    std::string tver(tokenver[n-1]);
    std::string tsim(tokensim[n-1]);
    
    for (int i = 0; i < std::min(tver.size(),tsim.size()); i++)
    {
      sec_ver.push_back(tver.at(i)-48);
      sec_sim.push_back(tsim.at(i)-48);
      lines.push_back(line_ver);
    }
             
    strcpy(bufver,"");
    strcpy(bufsim,"");
  }

  compa_files::do_ana(sec_ver,sec_sim,chipn);

  return;
}

//
// Here is the code doing the comparison between the two data streams, for a given CIC chip
//

void compa_files::do_ana(std::vector<int> sec_ver,std::vector<int> sec_sim, int n)
{
    std::cout << std::endl;
    std::cout << "Analizing chip " << n << std::endl;
    std::cout << std::endl;
    
    
    std::vector<int> start_ver; // Contains the starting bit
    std::vector<int> start_sim; // of of L1 words sent, along with the L1ID
    start_sim.clear();
    start_ver.clear();
    
    int new_seq = 0;
    int L1_id   = 0;

    //
    // First loops are intended to detect the start of events in the differents streams
    // along with the corresponding L1 ID
    //

    // Start with the verilog information

    // Look for CIC L1 headers (15 bits at 1)

    for (unsigned int i = 0; i < sec_ver.size(); i++)
    {
        if (new_seq==0 && sec_ver[i]==1)
        {
            new_seq++;
        }
        else if (new_seq!=0 && sec_ver[i]==1)
        {
            new_seq++;
        }
        else if (new_seq!=0 && sec_ver[i]==0)
        {
            new_seq=0;
        }
        
        if (new_seq==15) // This is SOE signal for the L1 CIC word
        {
            L1_id   = 0;
            new_seq = 0;
            
            for (int j = i+11; j < i+20; j++) // Get the L1 ID
            {
                L1_id = L1_id+sec_ver[j]*pow(2,i+19-j);
            }
            
            start_ver.push_back(i-15); // Starting bit of the L1 word
            start_ver.push_back(L1_id);// Corresponding L1 ID
        }
    }
    
    // Now do the same for simu stream

    L1_id   = 0;
    new_seq = 0;
    
    for (int i = 0; i < sec_sim.size(); i++) // Same for the sim info
    {
        if (new_seq==0 && sec_sim[i]==1)
        {
            new_seq++;
        }
        else if (new_seq!=0 && sec_sim[i]==1)
        {
            new_seq++;
        }
        else if (new_seq!=0 && sec_sim[i]==0)
        {
            new_seq=0;
        }
        
        if (new_seq==15)
        {
            L1_id   = 0;
            new_seq = 0;
            
            for (int j = i+11; j < i+20; j++)
            {
                L1_id = L1_id+sec_sim[j]*pow(2,i+19-j);
            }
            
            start_sim.push_back(i-15);
            start_sim.push_back(L1_id);
        }
    }
    
    //
    // Second loop, this is the complete analysis
    // 

    int index = std::min(start_ver.size()/2,start_sim.size()/2);
    
    int nmiss=0;
    int ireal=0;
    
    int error_bit_ver[9];
    int error_bit_sim[9];
    
    int nclus_ver[2];
    int nclus_sim[2];
    
    int rk;
    
    std::vector<int> ver_strip;
    std::vector<int> ver_pixel;
    
    std::vector<int> sim_strip;
    std::vector<int> sim_pixel;
    
    int npass = 0;
    
    for (int i = 0; i < index-1; i++) // Loop over events stored
    {
      // Logically there is currently no FIFO limitation in the simulation
      // so all sim L1 words are transmitted
      //
      // In the verilog the event might be skipped if the FIFO is full
      // we take care of that here

      while (start_sim[2*ireal+1] != start_ver[2*i+1])
      {
	++nmiss;
	std::cout << "!!! FIFO full for L1 ID " << start_sim[2*ireal+1] << std::endl;
	++ireal;
	if (ireal>=index-1) break;
      }
        
      if (ireal>=index-1) break;
        
      ++npass;

      std::cout << "Compare data content for L1 ID (VERILOG/SIMU) " << start_ver[2*i+1] << " / " << start_sim[2*ireal+1] << std::endl;
      std::cout << "-> Starting at BX " << start_ver[2*i]/8 << " / " << start_sim[2*ireal]/8 << std::endl;
        
      ver_pixel.clear();
      ver_strip.clear();
      sim_pixel.clear();
      sim_strip.clear();
      
      for (int j = 0; j < 9; j++)
      {
	error_bit_ver[j]=0;
	error_bit_sim[j]=0;
      }
      for (int j = 0; j < 2; j++)
      {
	nclus_ver[j]=0;
	nclus_sim[j]=0;
      }
        
      int nerr = 0;
        
      if (i<index-2) // Read the verilog data content
      {
	for (int j = 1; j < start_ver[2*i+2]-start_ver[2*i] ; j++)
	{
	  rk=j-1;
	  if (rk>=16 && rk<25) // Error info
	  {
	    error_bit_ver[rk%16] = sec_ver[start_ver[2*i]+j];
	    if (error_bit_ver[rk%16]) ++nerr;
	  }
                
	  if (rk==25 && nerr) // There is an error flag
	  {
	    for (int k = 0; k < 8; k++) // Bit 8 is for the CIC
	    {
	      if (error_bit_ver[k]) std::cout << "VER Chip " << k << " not received" << std::endl;
	    }
	  }
                
	  if (rk>=34 && rk<41) nclus_ver[0] += sec_ver[start_ver[2*i]+j]*pow(2,40-rk);
	  //	  if (rk==41 && nerr) std::cout << "N strip clusters " << nclus_ver[0] << std::endl;
	  if (rk==41) std::cout << "N strip clusters ver " << nclus_ver[0] << std::endl;
                
	  if (nclus_ver[0]!=0 && rk==41)
	  {
	    for (int k = 0; k < nclus_ver[0]; k++)
	    {
	      int chip = 0;
	      int pos  = 0;
	      int width= 0;
	      int mip  = 0;
              
	      int rkclus = start_ver[2*i]+42+7*m_isMPA+(13+m_isMPA)*k;
	      
	      for (int kk = rkclus;    kk < rkclus+3;  kk++) chip  += sec_ver[kk]*pow(2,rkclus+2-kk);
	      for (int kk = rkclus+3;  kk < rkclus+10; kk++) pos   += sec_ver[kk]*pow(2,rkclus+9-kk);
	      for (int kk = rkclus+10; kk < rkclus+13; kk++) width += sec_ver[kk]*pow(2,rkclus+12-kk);
	      if (m_isMPA) mip = sec_ver[rkclus+13];
              
	      ver_strip.push_back(chip);
	      ver_strip.push_back(pos);
	      ver_strip.push_back(width);
	      if (m_isMPA) ver_strip.push_back(mip);
              
	      //              std::cout << "Strip cluster " << k << " : " << chip << "/"
	      //                 << pos << "/" << width << "/" << mip << std::endl;
	    }
	  }
                
	  if (rk>=41 && rk<48 && m_isMPA) nclus_ver[1] += sec_ver[start_ver[2*i]+j]*pow(2,47-rk);
	  if (rk==48 && nerr && m_isMPA) std::cout << "N pixel clusters " << nclus_ver[1] << std::endl;
                
	  if (nclus_ver[1]!=0 && rk==48)
	  {
	    for (int k = 0; k < nclus_ver[1]; k++)
	    {
	      int chip = 0;
	      int pos  = 0;
	      int width= 0;
	      int z    = 0;
              
	      int rkclus = start_ver[2*i]+42+7*m_isMPA+(13+m_isMPA)*nclus_ver[0]+17*k;
              
	      for (int kk = rkclus;    kk < rkclus+3;  kk++) chip  += sec_ver[kk]*pow(2,rkclus+2-kk);
	      for (int kk = rkclus+3;  kk < rkclus+10; kk++) pos   += sec_ver[kk]*pow(2,rkclus+9-kk);
	      for (int kk = rkclus+10; kk < rkclus+13; kk++) width += sec_ver[kk]*pow(2,rkclus+12-kk);
	      for (int kk = rkclus+13; kk < rkclus+17; kk++) z     += sec_ver[kk]*pow(2,rkclus+16-kk);
              
	      ver_pixel.push_back(chip);
	      ver_pixel.push_back(pos);
	      ver_pixel.push_back(width);
	      ver_pixel.push_back(z);
              
	      //        std::cout << "Pixel cluster " << k << " : " << chip << "/"
	      //        << pos << "/" << width << "/" << z << std::endl;
	    }
	  }
	}
            
	nerr=0;
        
	for (int j = 1; j < start_sim[2*ireal+2]-start_sim[2*ireal] ; j++)
	{
	  rk=j-1;
	  if (rk>=16 && rk<25)
	  {
	    error_bit_sim[rk%16] = sec_sim[start_sim[2*ireal]+j];
	    if (error_bit_sim[rk%16]) ++nerr;
	  }
                
	  if (rk==25 && nerr)
	  {
	    for (int k = 0; k < 8; k++)
	      {
		if (error_bit_sim[k]) std::cout << "SIM Chip " << k << " not received" << std::endl;
	      }
	  }
                
	  if (rk>=34 && rk<41) nclus_sim[0] += sec_sim[start_sim[2*ireal]+j]*pow(2,40-rk);
	  //	  if (rk==41 && nerr) std::cout << "N strip clusters sim " << nclus_sim[0] << std::endl;
	  if (rk==41) std::cout << "N strip clusters sim " << nclus_sim[0] << std::endl;
                
	  if (nclus_sim[0]!=0 && rk==41)
	  {
	    for (int k = 0; k < nclus_sim[0]; k++)
	    {
	      int chip = 0;
	      int pos  = 0;
	      int width= 0;
	      int mip  = 0;
              
	      int rkclus = start_sim[2*ireal]+42+7*m_isMPA+(13+m_isMPA)*k;
              
	      for (int kk = rkclus;    kk < rkclus+3;  kk++) chip  += sec_sim[kk]*pow(2,rkclus+2-kk);
	      for (int kk = rkclus+3;  kk < rkclus+10; kk++) pos   += sec_sim[kk]*pow(2,rkclus+9-kk);
	      for (int kk = rkclus+10; kk < rkclus+13; kk++) width += sec_sim[kk]*pow(2,rkclus+12-kk);
	      if (m_isMPA) mip = sec_sim[rkclus+13];
              
	      sim_strip.push_back(chip);
	      sim_strip.push_back(pos);
	      sim_strip.push_back(width);
	      if (m_isMPA) sim_strip.push_back(mip);
              
	      //          std::cout << "Strip cluster " << k << " : " << chip << "/"
	      //         << pos << "/" << width << "/" << mip << std::endl;
	    }
	  }
                
                
	  if (rk>=41 && rk<48 && m_isMPA) nclus_sim[1] += sec_sim[start_sim[2*ireal]+j]*pow(2,47-rk);
	  if (rk==48 && nerr && m_isMPA) std::cout << "N pixel clusters " << nclus_sim[1] << std::endl;
          
	  if (nclus_sim[1]!=0 && rk==48)
	  {
	    for (int k = 0; k < nclus_sim[1]; k++)
	    {
	      int chip = 0;
	      int pos  = 0;
	      int width= 0;
	      int z    = 0;
              
	      int rkclus = start_sim[2*ireal]+42+7*m_isMPA+(13+m_isMPA)*nclus_sim[0]+17*k;
	      
	      for (int kk = rkclus;    kk < rkclus+3;  kk++) chip  += sec_sim[kk]*pow(2,rkclus+2-kk);
	      for (int kk = rkclus+3;  kk < rkclus+10; kk++) pos   += sec_sim[kk]*pow(2,rkclus+9-kk);
	      for (int kk = rkclus+10; kk < rkclus+13; kk++) width += sec_sim[kk]*pow(2,rkclus+12-kk);
	      for (int kk = rkclus+13; kk < rkclus+17; kk++) z     += sec_sim[kk]*pow(2,rkclus+16-kk);
              
	      sim_pixel.push_back(chip);
	      sim_pixel.push_back(pos);
	      sim_pixel.push_back(width);
	      sim_pixel.push_back(z);
              
	      //                        std::cout << "Pixel cluster " << k << " : " << chip << "/"
	      //                        << pos << "/" << width << "/" << z << std::endl;
	    }
	  }
	}
            
	if (nclus_ver[0]>0)
	{
	  for (unsigned int j = 0; j < ver_strip.size()/(3+m_isMPA) ; j++)
	  {
	    int chip = ver_strip.at((3+m_isMPA)*j);
            
	    if (error_bit_sim[chip]) continue;
            
	    int pos  = ver_strip.at((3+m_isMPA)*j+1);
	    int width= ver_strip.at((3+m_isMPA)*j+2);
	    int mip  = 0;

	    if (m_isMPA) mip = ver_strip.at((3+m_isMPA)*j+3);
            
	    bool found = false;
            
	    for (unsigned int jj = 0; jj < sim_strip.size()/(3+m_isMPA) ; jj++)
	    {
	      if (found) break;
              
	      if (chip != sim_strip.at((3+m_isMPA)*jj))   continue;
	      if (pos  != sim_strip.at((3+m_isMPA)*jj+1)) continue;
	      if (width!= sim_strip.at((3+m_isMPA)*jj+2)) continue;
	      if (m_isMPA && mip  != sim_strip.at((3+m_isMPA)*jj+3)) continue;
              
	      found = true;
	    }
            
	    if (!found)
	      std::cout << "ERROR!! -> Strip cluster " << j << " : " << chip << "/"
			<< pos << "/" << width << "/" << mip << " is in the verilog but not in the simu" << std::endl;
	  }
	}
        
	if (nclus_ver[1]>0)
	{
	  for (unsigned int j = 0; j < ver_pixel.size()/4 ; j++)
	  {
	    int chip = ver_pixel.at(4*j);
            
	    if (error_bit_sim[chip]) continue;
                    
	    int pos  = ver_pixel.at(4*j+1);
	    int width= ver_pixel.at(4*j+2);
	    int z    = ver_pixel.at(4*j+3);
            
	    bool found = false;
            
	    for (unsigned int jj = 0; jj < sim_pixel.size()/4 ; jj++)
	    {
	      if (found) break;
              
	      if (chip != sim_pixel.at(4*jj))   continue;
	      if (pos  != sim_pixel.at(4*jj+1)) continue;
	      if (width!= sim_pixel.at(4*jj+2)) continue;
	      if (z    != sim_pixel.at(4*jj+3)) continue;
              
	      found = true;
	    }
                    
	    if (!found)
	      std::cout << "ERROR!! -> Pixel cluster " << j << " : " << chip << "/"
                        << pos << "/" << width << "/" << z << " is in the verilog but not in the simu" << std::endl;
	  }
	}
            
	if (nclus_sim[0]>0)
	{
	  for (unsigned int j = 0; j < sim_strip.size()/(3+m_isMPA) ; j++)
	  {
	    int chip = sim_strip.at((3+m_isMPA)*j);
            
	    if (error_bit_ver[chip]) continue;
            
	    int pos  = sim_strip.at((3+m_isMPA)*j+1);
	    int width= sim_strip.at((3+m_isMPA)*j+2);
	    int mip  = 0;
	    
	    if (m_isMPA) mip  = sim_strip.at((3+m_isMPA)*j+3);
            
	    bool found = false;
	    
	    for (unsigned int jj = 0; jj < ver_strip.size()/(3+m_isMPA) ; jj++)
	    {
	      if (found) break;
              
	      if (chip != ver_strip.at((3+m_isMPA)*jj))   continue;
	      if (pos  != ver_strip.at((3+m_isMPA)*jj+1)) continue;
	      if (width!= ver_strip.at((3+m_isMPA)*jj+2)) continue;
	      if (m_isMPA && mip  != ver_strip.at((3+m_isMPA)*jj+3)) continue;
              
	      found = true;
	    }
            
	    if (!found)
	      std::cout << "ERROR!! -> Strip cluster " << j << " : " << chip << "/"
                        << pos << "/" << width << "/" << mip << " is in the simu but not in the verilog" << std::endl;
	  }
	}
            
	if (nclus_sim[1]>0)
	{
	  for (unsigned int j = 0; j < sim_pixel.size()/4 ; j++)
	  {
	    int chip = sim_pixel.at(4*j);
            
	    if (error_bit_ver[chip]) continue;
            
	    int pos  = sim_pixel.at(4*j+1);
	    int width= sim_pixel.at(4*j+2);
	    int z    = sim_pixel.at(4*j+3);
            
	    bool found = false;
            
	    for (unsigned int jj = 0; jj < ver_pixel.size()/4 ; jj++)
	    {
	      if (found) break;
              
	      if (chip != ver_pixel.at(4*jj))   continue;
	      if (pos  != ver_pixel.at(4*jj+1)) continue;
	      if (width!= ver_pixel.at(4*jj+2)) continue;
	      if (z    != ver_pixel.at(4*jj+3)) continue;
              
	      found = true;
	    }
            
	    if (!found)
	      std::cout << "ERROR!! -> Pixel cluster " << j << " : " << chip << "/"
                        << pos << "/" << width << "/" << z << " is in the simu but not in the verilog" << std::endl;
	  }
	}
            
	//for (int j = 1; j < start_ver[2*i+2]-start_ver[2*i] ; j++) std::cout << sec_ver[start_ver[2*i]+j];
	//std::cout << std::endl;
        
	//for (int j = 1; j < start_sim[2*ireal+2]-start_sim[2*ireal] ; j++) std::cout << sec_sim[start_sim[2*ireal]+j];
	//std::cout << std::endl;
	/*
	  for (int j = 1; j < start_ver[2*i+2]-start_ver[2*i] ; j++)
	  {
	  if (start_sim[2*ireal]+j>=start_sim[2*ireal+2]) continue;
	  if (sec_ver[start_ver[2*i]+j] != sec_sim[start_sim[2*ireal]+j])
	  std::cout << "Problem at BX " << (start_ver[2*i]+j)/8 << " / " << (start_sim[2*ireal]+j)/8 << std::endl;
	  }
	*/
      }
      ++ireal;
      /*
	if (sec_sim[i]!=sec_ver[i+5])
	{
	std::cout << lines[i] << "/" << sec_ver[i+5] << "/" << sec_sim[i] << std::endl;
	}
      */
    }
    std::cout << std::endl;
    
}

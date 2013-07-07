/*
 Copyright (c) 2012 Seyed Nematollah Ahmadyan
 All rights reserved.
 
 Developed by: 		Seyed Nematollah Ahmadyan [ahmadyan@gmail.com]
 Center of Reliability and High-Performance Computing,
 Coordinated Science Lab,
 Electrical and Computer Engineering Department,
 University of Illinois at Urbana-Champaign
 
 http://netfiles.uiuc.edu/ahmadya2/www
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal with the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
 Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimers in the documentation and/or other materials provided with the distribution.
 Neither the names of <Name of Development Group, Name of Institution>, nor the names of its contributors may be used to endorse or promote products derived from this Software without specific prior written permission.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
 */


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
#include <conio.h>   //for getch(), needed in wait_for_key()
#include <windows.h> //for Sleep()
void sleep(int i) { Sleep(i*1000); }
#endif


#define SLEEP_LGTH 2  // sleep time in seconds
#include "config.h"
#include "plot.h"
#include <sstream>
#include <cstdio>

Plotter::Plotter(string path, Configuration* config){
    gnuPlotPath = path ;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    string name;
    int name_enable ;
    config->getParameter("edu.uiuc.crhc.core.system.dimension", &dim);
    config->getParameter("edu.uiuc.crhc.core.system.x.min", &xmin);
    config->getParameter("edu.uiuc.crhc.core.system.x.max", &xmax);
    config->getParameter("edu.uiuc.crhc.core.system.y.min", &ymin);
    config->getParameter("edu.uiuc.crhc.core.system.y.max", &ymax);
    
    config->getParameter("edu.uiuc.crhc.core.options.nameEn", &name_enable);
    if(name_enable==1){
        config->getParameter("edu.uiuc.crhc.core.system.name", &name);
    }else{
        name = " " ;
    }
    
#ifdef _WIN32
    gnuplotPipe = _popen(gnuPlotPath.c_str(),"w");
#else
    gnuplotPipe = popen(gnuPlotPath.c_str(),"w");
#endif


    string buffer ;
    fflush(gnuplotPipe);
    cout << dim << " " << xmin << " " << xmax << endl ;
    emptyPlot(name, xmin, xmax, ymin, ymax);
    closed=false;
        
}

Plotter::~Plotter(){
    if(!closed) close();
}

void Plotter::close(){
    string buffer = "replot\n";
    fprintf(gnuplotPipe, buffer.c_str());
    fflush(gnuplotPipe);
    
    waitForKey();
    
#ifdef _WIN32
    _pclose(gnuplotPipe);
#else
    pclose(gnuplotPipe);
#endif


    closed=true;
}

void Plotter::emptyPlot(string title, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int angleX, int angleY, string xlabel, string ylabel, string zlabel){
    stringstream cmdstr;
    cmdstr << "splot [" << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "][" << zmin << ":" << zmax <<"] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << title << "\"  \n";
    cmdstr << "set xlabel \"$" << xlabel << "$\" \n";
    cmdstr << "set ylabel \"$" << ylabel << "$\" \n";
    cmdstr << "set zlabel \"$"<< zlabel<< "$\" \n";
    
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::emptyPlot(string title, double xmin, double xmax, double ymin, double ymax ){
    stringstream cmdstr;
    cmdstr << "plot [" << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << title << "\"  \n";
    
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}


void Plotter::drawLine(double iFromX, double iFromY, double iFromZ, double iToX, double iToY, double iToZ){
    stringstream cmdstr;
    cmdstr << " set arrow from " << iFromX << "," << iFromY << "," << iFromZ <<  "     to     " << iToX << "," << iToY  << "," << iToZ << "        nohead  lc rgb \"blue\" lw 2 \n" ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::drawLine(double iFromX, double iFromY, double iToX, double iToY){
    stringstream cmdstr;
    cmdstr << " set arrow from " << iFromX << "," << iFromY << " to " << iToX << "," << iToY << " nohead  lc rgb \"blue\" lw 2 \n" ;
    cout << cmdstr.str() << endl ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}


void Plotter::drawArray(vector< vector<double> > trace, int index1, int index2){
    if(trace.size()>0){
            //scaling
        
        double xmin=9999, xmax=-99999, ymin=9999, ymax=-9999;
        for(int i=0;i<trace.size();i++){
            if(trace.at(i).at(index1) > xmax) xmax=trace.at(i).at(index1);
            if(trace.at(i).at(index1) < xmin) xmin=trace.at(i).at(index1);
            if(trace.at(i).at(index2) > ymax) ymax=trace.at(i).at(index2);
            if(trace.at(i).at(index2) < ymin) ymin=trace.at(i).at(index2);

        }
        stringstream cmdstr;
        cmdstr << "plot [ " << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "] 0 with linespoints lt \"white\" pt 0.01" ;
        fprintf(gnuplotPipe, cmdstr.str().c_str());
        fflush(gnuplotPipe);

    for(int i=0;i<trace.size()-1;i++){
        vector<double> point1 = trace[i] ;
        vector<double> point2 = trace[i+1] ;
        drawLine( point1[index1], point1[index2], point2[index1], point2[index2] ) ;
    }
        
        
        stringstream cmdstr2;
        cmdstr2 << " replot \n" ;
        fprintf(gnuplotPipe, cmdstr2.str().c_str());
        fflush(gnuplotPipe);
    }
    else{ cout << "empty trace" << endl ;
    }
    
    
}

void Plotter::drawTrace(Trace* trace, int index1, int index2){
    drawArray(trace->getData(), index1, index2);
}


void Plotter::drawTrace(Trace* trace, double t_max, int index){
    
    double min=9999, max=-99999;
    for(int i=0;i<trace->getSampleSize();i++){
        if(trace->getSample(i)[index] > max) max=trace->getSample(i)[index];
        if(trace->getSample(i)[index] < min) min=trace->getSample(i)[index];
    }
    stringstream cmdstr;
    cmdstr << "plot [ " << 0 << ":" << t_max << "][" << min << ":" << max << "] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << " " << "\"  \n";
    //cmdstr << "set ylabel " << "y(" << index << " )"  << " \n";
    //cmdstr << "set xlabel " << "time" << " \n";
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
    
    double time_step = t_max / trace->getSampleSize();
    for(int i=0;i<trace->getData().size()-1;i++){
        vector<double> point1 = trace->getSample(i) ;
        vector<double> point2 = trace->getSample(i+1) ;
        drawLine( i*time_step, point1[index], (i+1)*time_step, point2[index] ) ;
    }
    stringstream cmdstr2;
    cmdstr2 << " replot \n" ;
    fprintf(gnuplotPipe, cmdstr2.str().c_str());
    fflush(gnuplotPipe);
}


void Plotter::waitForKey(){
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;
    
    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;
    
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

void Plotter::saveToPdf(string path){
    stringstream cmdstr;
    cmdstr << " set term post \n" ;
    cmdstr << " set output \"" << path << "\"\n" ;
    cmdstr << " replot \n" ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::execute(string str){
    fprintf(gnuplotPipe, str.c_str());
    fflush(gnuplotPipe);
    cout << str << endl ;
    string buffer = "replot\n";
    fprintf(gnuplotPipe, buffer.c_str());
    fflush(gnuplotPipe);
    
}

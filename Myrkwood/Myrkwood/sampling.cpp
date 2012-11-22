/*
 * sampling.cpp
 *
 *  Created on: Sep 19, 2012
 *      Author: adel
 */

#include "sampling.h"

sampling::sampling(Configuration* config) {

}

sampling::~sampling() {
}

System::Point* sampling::sample(){
	System::Point* p = new System::Point(2) ;
	return p;
}

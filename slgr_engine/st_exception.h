//
// Copyright 2016 MIT Lincoln Laboratory, Massachusetts Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use these files except in compliance with
// the License.
//
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.
//

//
// Exception for SLGR tools
//

#pragma once

class  ST_exception : public exception {
	string msg;
public:
 	ST_exception (const char *msg) {ST_exception::msg = msg;};
 	ST_exception (string &in_msg) {ST_exception::msg = in_msg;};
 	ST_exception (string in_msg) {ST_exception::msg = in_msg;};
 	ST_exception (const ST_exception &e) {msg = e.msg;};
 	string get_error_message (void) {return msg;};
	virtual ~ST_exception (void) throw() {};
	virtual const char* what() const throw() {
		return msg.c_str();
   };
};


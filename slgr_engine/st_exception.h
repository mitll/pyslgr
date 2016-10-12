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


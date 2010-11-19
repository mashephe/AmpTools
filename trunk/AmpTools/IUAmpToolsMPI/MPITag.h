#if !defined(MPITAG)
#define MPITAG


class MPITag {

 public:

  enum { kIntSend = 1,
	 kDoubleSend,
	 kCharSend,
	 kDataSend,
	 kAcknowledge,
	 kMaxTags };

};

#endif

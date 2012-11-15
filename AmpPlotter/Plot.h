#if !(defined PLOT)
#define PLOT

#include <string>

class Plot
{
	
public:
	
	Plot( std::string title, std::string dir, std::string tag ) :
	m_title( title ),
	m_directory( dir ),
	m_tag( tag ) {}
	
	std::string tag() const { return m_tag; }
	std::string directory() const { return m_directory; }
	std::string title() const { return m_title; }
	
private:
	
	std::string m_title;
	std::string m_directory;
	std::string m_tag;

};

#endif

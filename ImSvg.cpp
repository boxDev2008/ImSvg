#pragma once

#include "ImSvg.h"

#include "imgui_internal.h"
#include <fstream>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <functional>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ImSvg
{

struct XmlAttr { std::string name, value; };
struct XmlNode
{
    std::string tag;
    std::vector<XmlAttr> attrs;
    std::vector<XmlNode> children;
    std::string text;
    const std::string* Attr(const std::string& n) const
    {
        for (auto& a : attrs) if (a.name == n) return &a.value;
        return nullptr;
    }
};

struct SVGTransform
{
    float m[6] = {1,0,0,1,0,0};
    ImVec2 Apply(ImVec2 p) const
    {
        return { m[0]*p.x + m[2]*p.y + m[4], m[1]*p.x + m[3]*p.y + m[5] };
    }
    static SVGTransform Mul(const SVGTransform& A, const SVGTransform& B)
    {
        SVGTransform C;
        C.m[0]=A.m[0]*B.m[0]+A.m[2]*B.m[1]; C.m[1]=A.m[1]*B.m[0]+A.m[3]*B.m[1];
        C.m[2]=A.m[0]*B.m[2]+A.m[2]*B.m[3]; C.m[3]=A.m[1]*B.m[2]+A.m[3]*B.m[3];
        C.m[4]=A.m[0]*B.m[4]+A.m[2]*B.m[5]+A.m[4];
        C.m[5]=A.m[1]*B.m[4]+A.m[3]*B.m[5]+A.m[5];
        return C;
    }
};

enum class PathCmd { MoveTo, LineTo, CubicBez, QuadBez, Arc, Close };
struct PathSegment
{
    PathCmd cmd = PathCmd::MoveTo;
    ImVec2 p[3] = {};
    float rx=0, ry=0, xRot=0;
    bool largeArc=false, sweep=false;
};

namespace detail
{

inline void SkipWS(const char*& p)
{
    while (*p && (*p==' '||*p=='\t'||*p=='\n'||*p=='\r'||*p==',')) ++p;
}
inline void SkipWSNoComma(const char*& p)
{
    while (*p && (*p==' '||*p=='\t'||*p=='\n'||*p=='\r')) ++p;
}
inline float ParseFloat(const char*& p)
{
    SkipWS(p); char* e; float v=strtof(p,&e); p=e; return v;
}
inline bool ParseFlag(const char*& p)
{
    SkipWSNoComma(p);
    if (*p==',') ++p;
    SkipWSNoComma(p);
    bool v=(*p=='1'); ++p; return v;
}
inline ImVec2 ParseVec2(const char*& p)
{
    float x=ParseFloat(p),y=ParseFloat(p); return {x,y};
}

static const std::unordered_map<std::string,uint32_t> kColors =
{
    {"black",0x000000},{"white",0xffffff},{"red",0xff0000},{"green",0x008000},
    {"blue",0x0000ff},{"yellow",0xffff00},{"cyan",0x00ffff},{"magenta",0xff00ff},
    {"orange",0xffa500},{"purple",0x800080},{"gray",0x808080},{"grey",0x808080},
    {"silver",0xc0c0c0},{"maroon",0x800000},{"navy",0x000080},{"olive",0x808000},
    {"teal",0x008080},{"lime",0x00ff00},{"aqua",0x00ffff},{"fuchsia",0xff00ff},
    {"pink",0xffc0cb},{"brown",0xa52a2a},{"darkred",0x8b0000},{"transparent",0},
    {"lightblue",0xadd8e6},{"lightgreen",0x90ee90},{"lightyellow",0xffffe0},
    {"darkblue",0x00008b},{"darkgreen",0x006400},{"darkorange",0xff8c00},
    {"coral",0xff7f50},{"tomato",0xff6347},{"gold",0xffd700},{"khaki",0xf0e68c},
    {"crimson",0xdc143c},{"indigo",0x4b0082},{"violet",0xee82ee},{"beige",0xf5f5dc},
};

inline float HexV(char c)
{
    if(c>='0'&&c<='9')return(float)(c-'0');
    if(c>='a'&&c<='f')return(float)(c-'a'+10);
    if(c>='A'&&c<='F')return(float)(c-'A'+10);
    return 0;
}
inline float HexPair(char hi,char lo)
{
    return(HexV(hi)*16+HexV(lo))/255.f;
}

inline SVGColor ParseColor(const std::string& s)
{
    SVGColor c;
    if (s.empty()||s=="none"){c.none=true;return c;}
    if (s=="transparent"){c.none=true;return c;}
    if (s[0]=='#')
    {
        if (s.size()==7){c.r=HexPair(s[1],s[2]);c.g=HexPair(s[3],s[4]);c.b=HexPair(s[5],s[6]);}
        else if(s.size()==4){c.r=HexPair(s[1],s[1]);c.g=HexPair(s[2],s[2]);c.b=HexPair(s[3],s[3]);}
        else if(s.size()==9){c.r=HexPair(s[1],s[2]);c.g=HexPair(s[3],s[4]);c.b=HexPair(s[5],s[6]);c.a=HexPair(s[7],s[8]);}
        return c;
    }
    if (s.size()>5&&s.substr(0,5)=="rgba(")
    {
        const char* p=s.c_str()+5;
        c.r=ParseFloat(p)/255.f;c.g=ParseFloat(p)/255.f;c.b=ParseFloat(p)/255.f;c.a=ParseFloat(p);return c;
    }
    if (s.size()>4&&s.substr(0,4)=="rgb(")
    {
        const char* p=s.c_str()+4;
        c.r=ParseFloat(p)/255.f;c.g=ParseFloat(p)/255.f;c.b=ParseFloat(p)/255.f;return c;
    }
    auto it=kColors.find(s);
    if(it!=kColors.end())
    {
        uint32_t v=it->second;
        c.r=((v>>16)&0xff)/255.f;c.g=((v>>8)&0xff)/255.f;c.b=(v&0xff)/255.f;
    }
    return c;
}
inline float ParseLength(const std::string& s)
{
    const char* p=s.c_str();return ParseFloat(p);
}

inline void ApplyProp(const std::string& k,const std::string& v,SVGStyle& s)
{
    if(k=="fill")               s.fill         =ParseColor(v);
    else if(k=="stroke")        s.stroke       =ParseColor(v);
    else if(k=="stroke-width") {const char*p=v.c_str();s.strokeWidth =ParseFloat(p);}
    else if(k=="fill-opacity") {const char*p=v.c_str();s.fillOpacity =ParseFloat(p);}
    else if(k=="stroke-opacity"){const char*p=v.c_str();s.strokeOpacity=ParseFloat(p);}
    else if(k=="opacity")      {const char*p=v.c_str();s.opacity     =ParseFloat(p);}
}
inline void ParseStyleStr(const std::string& css,SVGStyle& s)
{
    std::istringstream ss(css);std::string tok;
    while(std::getline(ss,tok,';'))
    {
        auto c=tok.find(':');if(c==std::string::npos)continue;
        std::string k=tok.substr(0,c),v=tok.substr(c+1);
        auto trim=[](std::string& t)
        {
            while(!t.empty()&&(t.front()==' '||t.front()=='\t'))t.erase(t.begin());
            while(!t.empty()&&(t.back()==' '||t.back()=='\t'))t.pop_back();
        };
        trim(k);trim(v);ApplyProp(k,v,s);
    }
}
inline SVGStyle ParseStyle(const XmlNode& node,SVGStyle inh)
{
    SVGStyle s=inh;
    if(auto*v=node.Attr("fill"))          ApplyProp("fill",*v,s);
    if(auto*v=node.Attr("stroke"))        ApplyProp("stroke",*v,s);
    if(auto*v=node.Attr("stroke-width"))  ApplyProp("stroke-width",*v,s);
    if(auto*v=node.Attr("fill-opacity"))  ApplyProp("fill-opacity",*v,s);
    if(auto*v=node.Attr("stroke-opacity"))ApplyProp("stroke-opacity",*v,s);
    if(auto*v=node.Attr("opacity"))       ApplyProp("opacity",*v,s);
    if(auto*v=node.Attr("style"))         ParseStyleStr(*v,s);
    return s;
}

inline SVGTransform ParseTx(const std::string& s)
{
    SVGTransform t; const char* p=s.c_str();
    while(*p)
    {
        while(*p&&(*p==' '||*p=='\t'||*p=='\n'||*p=='\r'))++p;
        if(!*p)break;
        const char* start=p;while(*p&&*p!='(')++p;
        std::string fn(start,p);
        while(!fn.empty()&&fn.back()==' ')fn.pop_back();
        if(*p=='(')++p;
        SVGTransform cur;
        if(fn=="translate")
        {
            float tx=ParseFloat(p),ty=0;
            SkipWS(p);if(*p&&*p!=')')ty=ParseFloat(p);
            cur.m[4]=tx;cur.m[5]=ty;
        }
        else if(fn=="scale")
        {
            float sx=ParseFloat(p),sy=sx;
            SkipWS(p);if(*p&&*p!=')')sy=ParseFloat(p);
            cur.m[0]=sx;cur.m[3]=sy;
        }
        else if(fn=="rotate")
        {
            float a=ParseFloat(p)*(float)(M_PI/180.0),cx2=0,cy2=0;
            SkipWS(p);if(*p&&*p!=')'){cx2=ParseFloat(p);cy2=ParseFloat(p);}
            float co=cosf(a),si=sinf(a);
            cur.m[0]=co;cur.m[1]=si;cur.m[2]=-si;cur.m[3]=co;
            cur.m[4]=cx2-co*cx2+si*cy2;cur.m[5]=cy2-si*cx2-co*cy2;
        }
        else if(fn=="skewX"){float a=ParseFloat(p)*(float)(M_PI/180.0);cur.m[2]=tanf(a);}
        else if(fn=="skewY"){float a=ParseFloat(p)*(float)(M_PI/180.0);cur.m[1]=tanf(a);}
        else if(fn=="matrix")
        {
            cur.m[0]=ParseFloat(p);cur.m[1]=ParseFloat(p);
            cur.m[2]=ParseFloat(p);cur.m[3]=ParseFloat(p);
            cur.m[4]=ParseFloat(p);cur.m[5]=ParseFloat(p);
        }
        while(*p&&*p!=')')++p;if(*p==')')++p;
        t=SVGTransform::Mul(t,cur);
    }
    return t;
}

inline std::vector<PathSegment> ParsePath(const std::string& d)
{
    std::vector<PathSegment> segs;
    const char* p=d.c_str();
    char lastCmd=0;
    ImVec2 cur={0,0},startPt={0,0},lastCtrl={0,0};
    while(*p)
    {
        SkipWS(p);if(!*p)break;
        char cmd=*p;
        bool isLetter=(cmd>='A'&&cmd<='Z')||(cmd>='a'&&cmd<='z');
        if(isLetter){lastCmd=cmd;++p;}else cmd=lastCmd;
        SkipWS(p);
        bool isCurve=(cmd=='C'||cmd=='c'||cmd=='S'||cmd=='s'||
                      cmd=='Q'||cmd=='q'||cmd=='T'||cmd=='t');
        if(!isCurve) lastCtrl=cur;
        PathSegment s{};
        switch(cmd)
        {
        case 'M':case 'm':
        {
            ImVec2 pt=ParseVec2(p);
            if(cmd=='m')pt={cur.x+pt.x,cur.y+pt.y};
            cur=startPt=pt;
            s.cmd=PathCmd::MoveTo;s.p[0]=pt;segs.push_back(s);
            lastCmd=(cmd=='M')?'L':'l';break;
        }
        case 'L':case 'l':
        {
            ImVec2 pt=ParseVec2(p);if(cmd=='l')pt={cur.x+pt.x,cur.y+pt.y};
            s.cmd=PathCmd::LineTo;s.p[0]=cur;s.p[1]=pt;segs.push_back(s);cur=pt;break;
        }
        case 'H':case 'h':
        {
            float x=ParseFloat(p);if(cmd=='h')x+=cur.x;ImVec2 pt={x,cur.y};
            s.cmd=PathCmd::LineTo;s.p[0]=cur;s.p[1]=pt;segs.push_back(s);cur=pt;break;
        }
        case 'V':case 'v':
        {
            float y=ParseFloat(p);if(cmd=='v')y+=cur.y;ImVec2 pt={cur.x,y};
            s.cmd=PathCmd::LineTo;s.p[0]=cur;s.p[1]=pt;segs.push_back(s);cur=pt;break;
        }
        case 'C':case 'c':
        {
            ImVec2 cp1=ParseVec2(p),cp2=ParseVec2(p),ep=ParseVec2(p);
            if(cmd=='c'){cp1={cur.x+cp1.x,cur.y+cp1.y};cp2={cur.x+cp2.x,cur.y+cp2.y};ep={cur.x+ep.x,cur.y+ep.y};}
            s.cmd=PathCmd::CubicBez;s.p[0]=cur;s.p[1]=cp1;s.p[2]=ep;s.rx=cp2.x;s.ry=cp2.y;
            segs.push_back(s);lastCtrl=cp2;cur=ep;break;
        }
        case 'S':case 's':
        {
            ImVec2 cp2=ParseVec2(p),ep=ParseVec2(p);
            if(cmd=='s'){cp2={cur.x+cp2.x,cur.y+cp2.y};ep={cur.x+ep.x,cur.y+ep.y};}
            ImVec2 cp1={2*cur.x-lastCtrl.x,2*cur.y-lastCtrl.y};
            s.cmd=PathCmd::CubicBez;s.p[0]=cur;s.p[1]=cp1;s.p[2]=ep;s.rx=cp2.x;s.ry=cp2.y;
            segs.push_back(s);lastCtrl=cp2;cur=ep;break;
        }
        case 'Q':case 'q':
        {
            ImVec2 cp=ParseVec2(p),ep=ParseVec2(p);
            if(cmd=='q'){cp={cur.x+cp.x,cur.y+cp.y};ep={cur.x+ep.x,cur.y+ep.y};}
            s.cmd=PathCmd::QuadBez;s.p[0]=cur;s.p[1]=cp;s.p[2]=ep;
            segs.push_back(s);lastCtrl=cp;cur=ep;break;
        }
        case 'T':case 't':
        {
            ImVec2 ep=ParseVec2(p);if(cmd=='t')ep={cur.x+ep.x,cur.y+ep.y};
            ImVec2 cp={2*cur.x-lastCtrl.x,2*cur.y-lastCtrl.y};
            s.cmd=PathCmd::QuadBez;s.p[0]=cur;s.p[1]=cp;s.p[2]=ep;
            segs.push_back(s);lastCtrl=cp;cur=ep;break;
        }
        case 'A':case 'a':
        {
            s.cmd=PathCmd::Arc;
            s.rx=ParseFloat(p);s.ry=ParseFloat(p);s.xRot=ParseFloat(p);
            s.largeArc=ParseFlag(p);s.sweep=ParseFlag(p);
            ImVec2 ep=ParseVec2(p);if(cmd=='a')ep={cur.x+ep.x,cur.y+ep.y};
            s.p[0]=cur;s.p[2]=ep;segs.push_back(s);cur=ep;break;
        }
        case 'Z':case 'z':
        {
            s.cmd=PathCmd::Close;s.p[0]=cur;s.p[1]=startPt;segs.push_back(s);cur=startPt;break;
        }
        default:++p;break;
        }
    }
    return segs;
}

inline void XmlSkipWS(const char*& p)
{
    while(*p&&(*p==' '||*p=='\t'||*p=='\n'||*p=='\r'))++p;
}
inline std::string XmlReadVal(const char*& p)
{
    char q=*p++;std::string v;
    while(*p&&*p!=q)
    {
        if(*p=='&')
        {
            const char*e=p+1;while(*e&&*e!=';')++e;std::string en(p+1,e);
            if(en=="amp")v+='&';else if(en=="lt")v+='<';else if(en=="gt")v+='>';
            else if(en=="quot")v+='"';else if(en=="apos")v+='\'';
            p=(*e?e+1:e);
        }
        else v+=*p++;
    }
    if(*p)++p;return v;
}
inline std::string XmlReadName(const char*& p)
{
    std::string n;while(*p&&*p!=' '&&*p!='\t'&&*p!='\n'&&*p!='\r'&&*p!='='&&*p!='>'&&*p!='/')n+=*p++;return n;
}
bool ParseXmlNode(const char*& p, XmlNode& node);
inline void ParseXmlChildren(const char*& p, XmlNode& parent)
{
    while(*p)
    {
        XmlSkipWS(p);if(!*p)break;
        if(*p=='<')
        {
            if(*(p+1)=='/')break;
            if(*(p+1)=='!')
            {
                if(strncmp(p,"<!--",4)==0){const char*e=strstr(p,"-->");p=(e?e+3:p+4);continue;}
                while(*p&&*p!='>')++p;if(*p)++p;continue;
            }
            if(*(p+1)=='?'){while(*p&&*p!='>')++p;if(*p)++p;continue;}
            XmlNode child;if(ParseXmlNode(p,child))parent.children.push_back(std::move(child));
        }
        else{while(*p&&*p!='<')parent.text+=*p++;}
    }
}
inline bool ParseXmlNode(const char*& p, XmlNode& node)
{
    if(*p!='<')return false;
    ++p;XmlSkipWS(p);node.tag=XmlReadName(p);
    while(true)
    {
        XmlSkipWS(p);if(!*p||*p=='>'||(*p=='/'&&*(p+1)=='>'))break;
        std::string name=XmlReadName(p);XmlSkipWS(p);std::string value;
        if(*p=='='){++p;XmlSkipWS(p);if(*p=='"'||*p=='\'')value=XmlReadVal(p);}
        if(!name.empty())node.attrs.push_back({name,value});
    }
    if(*p=='/'){++p;if(*p=='>')++p;return true;}
    if(*p=='>')++p;
    ParseXmlChildren(p,node);
    if(*p=='<'&&*(p+1)=='/'){while(*p&&*p!='>')++p;if(*p)++p;}
    return true;
}
inline XmlNode ParseXml(const std::string& xml)
{
    XmlNode root;const char* p=xml.c_str();
    while(*p)
    {
        XmlSkipWS(p);
        if(*p=='<'&&*(p+1)=='?'){while(*p&&*p!='>')++p;if(*p)++p;continue;}
        if(*p=='<'&&*(p+1)=='!'){while(*p&&*p!='>')++p;if(*p)++p;continue;}
        break;
    }
    ParseXmlNode(p,root);return root;
}

struct Bez4{ImVec2 p0,p1,p2,p3;};
inline void ArcToBeziers(ImVec2 from,ImVec2 to,float rx,float ry, float xRotDeg,bool largeArc,bool sweep,std::vector<Bez4>& out)
{
    if(rx==0||ry==0)return;
    float phi=xRotDeg*(float)(M_PI/180.0);
    float cosPhi=cosf(phi),sinPhi=sinf(phi);
    float dx=(from.x-to.x)*0.5f,dy=(from.y-to.y)*0.5f;
    float x1p= cosPhi*dx+sinPhi*dy;
    float y1p=-sinPhi*dx+cosPhi*dy;
    rx=fabsf(rx);ry=fabsf(ry);
    float rx2=rx*rx,ry2=ry*ry,x1p2=x1p*x1p,y1p2=y1p*y1p;
    float lam=sqrtf(x1p2/rx2+y1p2/ry2);
    if(lam>1.f){rx*=lam;ry*=lam;rx2=rx*rx;ry2=ry*ry;}
    float num=rx2*ry2-rx2*y1p2-ry2*x1p2;
    float den=rx2*y1p2+ry2*x1p2;
    float sq=(den==0.f)?0.f:sqrtf(fabsf(num/den));
    if(largeArc==sweep)sq=-sq;
    float cxp= sq*rx*y1p/ry;
    float cyp=-sq*ry*x1p/rx;
    float cx=cosPhi*cxp-sinPhi*cyp+(from.x+to.x)*0.5f;
    float cy=sinPhi*cxp+cosPhi*cyp+(from.y+to.y)*0.5f;
    auto ang=[](float ux,float uy,float vx,float vy)->float
    {
        float n=sqrtf((ux*ux+uy*uy)*(vx*vx+vy*vy));if(n==0)return 0.f;
        float c=fmaxf(-1.f,fminf(1.f,(ux*vx+uy*vy)/n));float a=acosf(c);
        return(ux*vy-uy*vx<0)?-a:a;
    };
    float theta1=ang(1,0,(x1p-cxp)/rx,(y1p-cyp)/ry);
    float dtheta=ang((x1p-cxp)/rx,(y1p-cyp)/ry,(-x1p-cxp)/rx,(-y1p-cyp)/ry);
    if(!sweep&&dtheta>0)dtheta-=2.f*(float)M_PI;
    if( sweep&&dtheta<0)dtheta+=2.f*(float)M_PI;
    int nSeg=(int)ceilf(fabsf(dtheta)/(float)(M_PI*0.5f));if(nSeg==0)return;
    float dt=dtheta/nSeg;
    float tdt2=tanf(dt*0.5f);
    float alphaC=(fabsf(dt)<1e-5f)?0.f:sinf(dt)*(sqrtf(4.f+3.f*tdt2*tdt2)-1.f)/3.f;
    auto ptE=[&](float t)->ImVec2{float co=cosf(t),si=sinf(t);return{cosPhi*rx*co-sinPhi*ry*si+cx,sinPhi*rx*co+cosPhi*ry*si+cy};};
    auto dpE=[&](float t)->ImVec2{float co=cosf(t),si=sinf(t);return{cosPhi*rx*(-si)-sinPhi*ry*co,sinPhi*rx*(-si)+cosPhi*ry*co};};
    float t=theta1;ImVec2 P1=ptE(t),dP1=dpE(t);
    for(int i=0;i<nSeg;++i)
    {
        t+=dt;ImVec2 P2=ptE(t),dP2=dpE(t);
        Bez4 b;b.p0=P1;b.p1={P1.x+alphaC*dP1.x,P1.y+alphaC*dP1.y};
        b.p2={P2.x-alphaC*dP2.x,P2.y-alphaC*dP2.y};b.p3=P2;out.push_back(b);
        P1=P2;dP1=dP2;
    }
}

static const int kCS=20, kQS=16, kAS=20;
inline void TessPath(const std::vector<PathSegment>& path, const SVGTransform& xf, std::vector<Contour>& out)
{
    Contour cur;
    auto flush=[&]{ if(!cur.pts.empty()){out.push_back(cur);cur.pts.clear();cur.closed=false;} };
    auto T=[&](ImVec2 p)->ImVec2{ return xf.Apply(p); };
    for(auto& seg:path)
    {
        switch(seg.cmd)
        {
        case PathCmd::MoveTo: flush(); cur.pts.push_back(T(seg.p[0])); break;
        case PathCmd::LineTo: cur.pts.push_back(T(seg.p[1])); break;
        case PathCmd::CubicBez:
        {
            ImVec2 p0=T(seg.p[0]),p1=T(seg.p[1]),p2=T({seg.rx,seg.ry}),p3=T(seg.p[2]);
            for(int i=1;i<=kCS;++i)
            {
                float t=(float)i/kCS,mt=1-t,mt2=mt*mt,t2=t*t,mt3=mt2*mt,t3=t2*t;
                cur.pts.push_back({mt3*p0.x+3*mt2*t*p1.x+3*mt*t2*p2.x+t3*p3.x,
                                   mt3*p0.y+3*mt2*t*p1.y+3*mt*t2*p2.y+t3*p3.y});
            }
            break;
        }
        case PathCmd::QuadBez:
        {
            ImVec2 p0=T(seg.p[0]),p1=T(seg.p[1]),p2=T(seg.p[2]);
            for(int i=1;i<=kQS;++i)
            {
                float t=(float)i/kQS,mt=1-t;
                cur.pts.push_back({mt*mt*p0.x+2*mt*t*p1.x+t*t*p2.x,
                                   mt*mt*p0.y+2*mt*t*p1.y+t*t*p2.y});
            }
            break;
        }
        case PathCmd::Arc:
        {
            std::vector<Bez4> bzs;
            ArcToBeziers(seg.p[0],seg.p[2],seg.rx,seg.ry,seg.xRot,seg.largeArc,seg.sweep,bzs);
            for(auto& bz:bzs)
            {
                ImVec2 p0=T(bz.p0),p1=T(bz.p1),p2=T(bz.p2),p3=T(bz.p3);
                for(int i=1;i<=kAS;++i)
                {
                    float t=(float)i/kAS,mt=1-t,mt2=mt*mt,t2=t*t,mt3=mt2*mt,t3=t2*t;
                    cur.pts.push_back({mt3*p0.x+3*mt2*t*p1.x+3*mt*t2*p2.x+t3*p3.x,
                                       mt3*p0.y+3*mt2*t*p1.y+3*mt*t2*p2.y+t3*p3.y});
                }
            }
            break;
        }
        case PathCmd::Close: cur.closed=true; flush(); break;
        }
    }
    flush();
}

inline float Cross2(ImVec2 O,ImVec2 A,ImVec2 B)
{
    return(A.x-O.x)*(B.y-O.y)-(A.y-O.y)*(B.x-O.x);
}
inline bool PtInTri(ImVec2 p,ImVec2 a,ImVec2 b,ImVec2 c)
{
    float d1=Cross2(p,a,b),d2=Cross2(p,b,c),d3=Cross2(p,c,a);
    return !((d1<0||d2<0||d3<0)&&(d1>0||d2>0||d3>0));
}
inline float PolyArea(const std::vector<ImVec2>& v)
{
    float a=0;int n=(int)v.size();
    for(int i=0;i<n;++i){int j=(i+1)%n;a+=v[i].x*v[j].y-v[j].x*v[i].y;}
    return a*0.5f;
}
inline std::vector<Triangle> EarClip(const std::vector<ImVec2>& poly)
{
    std::vector<Triangle> tris;
    int n=(int)poly.size();if(n<3)return tris;
    if(n==3){tris.push_back({0,1,2});return tris;}
    std::vector<int> idx(n);for(int i=0;i<n;++i)idx[i]=i;
    if(PolyArea(poly)<0)std::reverse(idx.begin(),idx.end());
    int limit=n*n+16;
    while((int)idx.size()>3&&limit-->0)
    {
        int m=(int)idx.size();bool found=false;
        for(int i=0;i<m;++i)
        {
            int prev=(i-1+m)%m,next=(i+1)%m;
            ImVec2 A=poly[idx[prev]],B=poly[idx[i]],C=poly[idx[next]];
            if(Cross2(A,B,C)<=1e-6f)continue;
            bool ear=true;
            for(int j=0;j<m&&ear;++j)
            {
                if(j==prev||j==i||j==next)continue;
                if(PtInTri(poly[idx[j]],A,B,C))ear=false;
            }
            if(ear){tris.push_back({idx[prev],idx[i],idx[next]});idx.erase(idx.begin()+i);found=true;break;}
        }
        if(!found)break;
    }
    if((int)idx.size()==3)tris.push_back({idx[0],idx[1],idx[2]});
    return tris;
}

inline void CollectShapes(const XmlNode& node,SVGDocument& doc, SVGStyle inh,SVGTransform parentXf)
{
    const std::string& tag=node.tag;
    if(tag=="defs"||tag=="metadata"||tag=="title"||tag=="desc"||tag=="symbol")return;
    SVGStyle style=ParseStyle(node,inh);
    SVGTransform xf=parentXf;
    if(auto*v=node.Attr("transform"))xf=SVGTransform::Mul(xf,ParseTx(*v));
    auto getF=[&](const char* n,float def=0)->float{if(auto*v=node.Attr(n))return ParseLength(*v);return def;};

    if(tag=="g"||tag=="svg"||tag=="a")
    {
        for(auto& c:node.children)CollectShapes(c,doc,style,xf);return;
    }

    RenderShape rs; rs.style=style;

    auto bakeContour=[&](Contour c)
    {
        std::vector<Triangle> tris;
        if(c.closed&&c.pts.size()>=3)tris=EarClip(c.pts);
        rs.contours.push_back(std::move(c));
        rs.fillTris.push_back(std::move(tris));
    };

    if(tag=="path")
    {
        if(auto*v=node.Attr("d"))
        {
            auto segs=ParsePath(*v);
            TessPath(segs,xf,rs.contours);
            for(auto& c:rs.contours)
            {
                std::vector<Triangle> tris;
                if(c.closed&&c.pts.size()>=3)tris=EarClip(c.pts);
                rs.fillTris.push_back(std::move(tris));
            }
        }
    }
    else if(tag=="rect")
    {
        float x=getF("x"),y=getF("y"),w=getF("width"),h=getF("height");
        float rx=getF("rx"),ry2=getF("ry");
        if(rx>0&&ry2==0)ry2=rx;if(ry2>0&&rx==0)rx=ry2;
        Contour c;c.closed=true;
        if(rx<=0.5f)
        {
            c.pts.push_back(xf.Apply({x,y}));c.pts.push_back(xf.Apply({x+w,y}));
            c.pts.push_back(xf.Apply({x+w,y+h}));c.pts.push_back(xf.Apply({x,y+h}));
        }
        else
        {
            int seg=8;
            auto corner=[&](float ox,float oy,float startA)
            {
                for(int i=0;i<=seg;++i){float a=startA+(float)M_PI*0.5f*i/seg;c.pts.push_back(xf.Apply({ox+rx*cosf(a),oy+ry2*sinf(a)}));}
            };
            corner(x+w-rx,y+ry2,  -(float)M_PI*0.5f);
            corner(x+rx,  y+ry2,  -(float)M_PI);
            corner(x+rx,  y+h-ry2,-(float)M_PI*1.5f);
            corner(x+w-rx,y+h-ry2, 0.f);
        }
        bakeContour(std::move(c));
    }
    else if(tag=="circle")
    {
        float cx2=getF("cx"),cy2=getF("cy"),r=getF("r");
        Contour c;c.closed=true;int N=64;
        for(int i=0;i<N;++i){float a=2.f*(float)M_PI*i/N;c.pts.push_back(xf.Apply({cx2+r*cosf(a),cy2+r*sinf(a)}));}
        bakeContour(std::move(c));
    }
    else if(tag=="ellipse")
    {
        float cx2=getF("cx"),cy2=getF("cy"),ex=getF("rx"),ey=getF("ry");
        Contour c;c.closed=true;int N=64;
        for(int i=0;i<N;++i){float a=2.f*(float)M_PI*i/N;c.pts.push_back(xf.Apply({cx2+ex*cosf(a),cy2+ey*sinf(a)}));}
        bakeContour(std::move(c));
    }
    else if(tag=="line")
    {
        rs.isLine=true;rs.lineA=xf.Apply({getF("x1"),getF("y1")});rs.lineB=xf.Apply({getF("x2"),getF("y2")});
    }
    else if(tag=="polyline"||tag=="polygon")
    {
        Contour c;c.closed=(tag=="polygon");
        if(auto*v=node.Attr("points")){const char*p=v->c_str();while(*p){SkipWS(p);if(!*p)break;c.pts.push_back(xf.Apply(ParseVec2(p)));}}
        bakeContour(std::move(c));
    }
    else {return;}

    doc.shapes.push_back(std::move(rs));
}

}

SVGDocument ParseSVG(const std::string& xml)
{
    SVGDocument doc;
    auto root=detail::ParseXml(xml);
    const XmlNode* svgNode=nullptr;
    if(root.tag=="svg")svgNode=&root;
    else for(auto& c:root.children)if(c.tag=="svg"){svgNode=&c;break;}
    if(!svgNode)return doc;
    auto getF=[&](const char* n,float def=0)->float{if(auto*v=svgNode->Attr(n))return detail::ParseLength(*v);return def;};
    doc.width=getF("width",100);doc.height=getF("height",100);
    if(auto*vb=svgNode->Attr("viewBox"))
    {
        const char*p=vb->c_str();
        doc.viewBox[0]=detail::ParseFloat(p);doc.viewBox[1]=detail::ParseFloat(p);
        doc.viewBox[2]=detail::ParseFloat(p);doc.viewBox[3]=detail::ParseFloat(p);
        doc.hasViewBox=true;
    }
    SVGStyle def;def.fill={0,0,0,1,false};def.stroke.none=true;
    SVGTransform identity;
    for(auto& c:svgNode->children)detail::CollectShapes(c,doc,def,identity);
    return doc;
}

SVGDocument LoadSVG(const char* path)
{
    std::ifstream f(path);
    std::string xml((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    f.close();
    return ParseSVG(xml);
}

void RenderSVG(ImDrawList* dl, const SVGDocument& doc, ImVec2 origin, float scale, float globalAlpha)
{
    if(!dl) return;
    float vbSX=1, vbSY=1, vbOX=0, vbOY=0;
    if(doc.hasViewBox && doc.viewBox[2]>0 && doc.viewBox[3]>0)
    {
        vbSX=doc.width/doc.viewBox[2];  vbSY=doc.height/doc.viewBox[3];
        vbOX=-doc.viewBox[0]*vbSX;      vbOY=-doc.viewBox[1]*vbSY;
    }
    auto S=[&](ImVec2 p)->ImVec2
    {
        return{ origin.x+(p.x*vbSX+vbOX)*scale,
                origin.y+(p.y*vbSY+vbOY)*scale };
    };
    float strokeScale=scale*(vbSX+vbSY)*0.5f;

    auto mkCol=[&](const SVGColor& c, float extraAlpha)->ImU32
    {
        float a=fmaxf(0.f, fminf(1.f, c.a * extraAlpha * globalAlpha));
        return IM_COL32((int)(c.r*255.f+0.5f),(int)(c.g*255.f+0.5f),
                        (int)(c.b*255.f+0.5f),(int)(a*255.f+0.5f));
    };

    for(auto& rs : doc.shapes)
    {
        const SVGStyle& st = rs.style;
        float effFillA   = st.opacity * st.fillOpacity;
        float effStrokeA = st.opacity * st.strokeOpacity;
        bool hasFill  = !st.fill.none   && st.fill.a   > 0 && effFillA   > 0;
        bool hasStroke= !st.stroke.none && st.stroke.a > 0 && effStrokeA > 0 && st.strokeWidth > 0;

        ImU32 fillCol   = hasFill   ? mkCol(st.fill,   effFillA)   : 0;
        ImU32 strokeCol = hasStroke ? mkCol(st.stroke, effStrokeA) : 0;
        float sw = st.strokeWidth * strokeScale;

        if(rs.isLine)
        {
            if(hasStroke) dl->AddLine(S(rs.lineA), S(rs.lineB), strokeCol, sw);
            continue;
        }

        if(hasFill)
        {
            for(int ci=0; ci<(int)rs.contours.size(); ++ci)
            {
                auto& pts = rs.contours[ci].pts;
                if(ci >= (int)rs.fillTris.size()) break;
                auto& tris = rs.fillTris[ci];
                if(tris.empty()) continue;

                int vtxCount = (int)pts.size();
                int idxCount = (int)tris.size() * 3;
                dl->PrimReserve(idxCount, vtxCount);

                ImDrawIdx vtxStart = (ImDrawIdx)dl->_VtxCurrentIdx;
                for(auto& p : pts)
                {
                    dl->_VtxWritePtr->pos = S(p);
                    dl->_VtxWritePtr->uv  = dl->_Data->TexUvWhitePixel;
                    dl->_VtxWritePtr->col = fillCol;
                    dl->_VtxWritePtr++;
                    dl->_VtxCurrentIdx++;
                }
                for(auto& tri : tris)
                {
                    dl->_IdxWritePtr[0] = (ImDrawIdx)(vtxStart + tri.a);
                    dl->_IdxWritePtr[1] = (ImDrawIdx)(vtxStart + tri.b);
                    dl->_IdxWritePtr[2] = (ImDrawIdx)(vtxStart + tri.c);
                    dl->_IdxWritePtr += 3;
                }
            }
            dl->AddDrawCmd();
        }

        if(hasStroke)
        {
            for(auto& c : rs.contours)
            {
                if(c.pts.empty()) continue;
                dl->PathClear();
                for(auto& p : c.pts) dl->PathLineToMergeDuplicate(S(p));
                if(c.closed) dl->PathStroke(strokeCol, ImDrawFlags_Closed, sw);
                else         dl->PathStroke(strokeCol, ImDrawFlags_None,   sw);
            }
            dl->AddDrawCmd();
        }
    }
}

void RenderSVGString(ImDrawList* dl,const std::string& svgXml, ImVec2 origin,float scale,float globalAlpha)
{
    RenderSVG(dl,ParseSVG(svgXml),origin,scale,globalAlpha);
}

float SvgFitScale(const SVGDocument& doc,float boxW,float boxH)
{
    if(doc.width<=0||doc.height<=0)return 1.f;
    return fminf(boxW/doc.width,boxH/doc.height);
}

}
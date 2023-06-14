from generate_data import parse_input, extract_data
from writing import write_data
from art import text2art

def main():

    grouped = parse_input()
    gene_list = extract_data(grouped)
    write_data(gene_list)
    return 0

if __name__ == "__main__":
    
    Art=text2art("ORFribo+",font='block',chr_ignore=True) 
    print(Art)

    main()

    Art=text2art("RUN FINISHED",font='block',chr_ignore=True) 
    print(Art)  

